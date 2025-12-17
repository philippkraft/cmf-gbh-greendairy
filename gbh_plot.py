import tables
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from gbh_model import load_soilmoisture
from contextlib import contextmanager


@contextmanager
def get_data(code):
    with tables.open_file(f'gbh_{code}.h5') as res_file:
        yield res_file.get_node('/data')


def get_range(code, data):
    p_names = data.colnames[6:-1]
    print(f'n_total =', len(data.col('RMSEcalib')))
    print(f'RMSEmin = {np.min(data.col("RMSEcalib")):0.4f}')
    threshold, take = get_threshold(code, data)    
    print(f'n = {np.sum(take)} for RMSE<{threshold:0.4f}')
    if np.any(take):
        pn_len = max(len(pn) for pn in p_names)
        for pn in p_names:
            pmin = np.percentile(data.col(pn)[take], 5)
            pmax = np.percentile(data.col(pn)[take], 95)
            print(f'{code} - {pn.ljust(pn_len)}: {pmin:0.3f} - {pmax:0.3f}')
    return threshold

def get_threshold(code, data):
    for threshold in np.arange(0.01,0.1,0.005):
        take = np.array(data.col('RMSEcalib')[:] < threshold)
        if np.sum(take) > 10:
            return threshold, take
    else:
        raise ValueError('No or < 10 data points with RMSE<0.1 found in model ' + code)

def get_tabdata(code, data, threshold):
    take = data.col('RMSEcalib')[:] < threshold
    n = np.sum(take)
    def make_row(valuetype):
        col = data.col(valuetype)[take]
        return {
            'code': code,
            'n': n,
            'threshold': threshold,
            'type': valuetype,
            'p5': np.percentile(col, 5),
            'mean' : np.mean(col),
            'p95' : np.percentile(col, 95)
         }
    return code, [
        make_row(cn)
        for cn in data.colnames[:-1]
    ]

def dotty_plot(code, data, threshold):
    """
    Plots the distribution of parameters in the behavioural runs
    """
    take = np.array(data.col('RMSEcalib')[:] < threshold * 2)
    fig = plt.figure(figsize=(8, 8), dpi=100)
    p_names = data.colnames[6:-1]
    for i, pn in enumerate(p_names):
        plt.subplot(((len(p_names)-1) // 3) + 1, 3, i + 1)
        params = data.col(pn)[take]
        plt.plot(params, -data.col('RMSEcalib')[take], '.', alpha=1/np.log10(np.sum(take)))
        # plt.title(pn.replace('_', ' '), loc='left', fontsize=10)
        try:
            gk = gaussian_kde(params)
            x = np.linspace(params.min(), params.max(), 1001)
            plt.ylabel('RMSEcalib')
            plt.xlabel(pn.replace('_', ' '))
            plt.twinx()
            plt.plot(x, gk(x), 'r:')
            plt.yticks([])

        except (np.linalg.linalg.LinAlgError, ValueError):
            pass
    plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.4)
    plt.suptitle(f'Model {code}, RMSEcalib<{threshold:0.3f}')
    plt.grid()

    return fig

def plot_equifinality(code, data, threshold):
    take = np.array(data.col('RMSEcalib')[:] < threshold)
    r_names = data.colnames[2:6]
    fig = plt.figure(figsize=(8, 4), dpi=100)
    plt.violinplot([data.col(rn)[take] for rn in r_names])
    plt.xticks(range(1, len(r_names)+1), r_names)        
    plt.title(f'Model {code}, RMSEcalib<{threshold:0.3f}')
    plt.grid()
    return fig

def plot_validation(code, data, threshold):
    take = np.array(data.col('RMSEcalib')[:] < threshold)
    fig = plt.figure(figsize=(8, 8), dpi=100)

    plt.scatter(data.col('RMSEcalib')[take], data.col('RMSEvalid')[take])
    plt.plot([0,.1],[0,.1], 'k:')
    plt.title(f'Model {code}, RMSEcalib<{threshold:0.3f}')
    plt.grid()
    return fig

def plot_moisture(code, data, threshold):
    plot_hours = 3*365*24
    fig = plt.figure(figsize=(8,4), dpi=100)
    sim_index = pd.date_range('2021-11-01','2025-08-01', freq='h')[:-1]
    take = np.array(data.col('RMSEcalib')[:] < threshold)
    for sim in data.col('simulation')[take, -plot_hours:]:
        plt.plot(sim_index[-plot_hours:], sim[:,-2], color='blue', alpha=0.1)
        plt.plot(sim_index[-plot_hours:], sim[:,-1], color='orange', alpha=0.1)
        plt.title(f'Model {code}, RMSEcalib<{threshold:0.3f}')
    
    for observation in load_soilmoisture('moisture.xlsx'):
        observation = observation[[c for c in observation.columns if 'mean' in c]]
        observation.plot(linestyle='--', ax=plt.gca())
    return fig

def plot_boundary(code, data, threshold):
    fig = plt.figure(figsize=(8,4), dpi=100)
    take = np.array(data.col('RMSEcalib')[:] < threshold)
    sim = np.cumsum(data.col('simulation')[take,-365*24:,:4].T/24, 1)
    sim_index = pd.date_range('2021-11-01','2024-01-01', freq='h')[:-1][-365*24:]    
    plt.plot(sim_index, sim[3], alpha=0.3, label='RO', c='r')
    plt.plot(sim_index, sim[0], alpha=0.3, label='P', c='b') 
    plt.plot(sim_index, sim[1], alpha=0.3, label='ET', c='g')
    plt.plot(sim_index, sim[2], alpha=0.3, label='GW', c='c')
    plt.ylabel('mm')
    plt.title(f'Model {code}, RMSEcalib<{threshold:0.3f}')
    return fig

def result(code, data, threshold):
    """
    Creates a DataFrame with summary statistics of the behavioural simulations
    
    """
    take = np.array(data.col('RMSEcalib')[:] < threshold)
    sim = data.col('simulation')[take]
    sim_index = pd.date_range('2020-08-26 12:00:00', periods=sim.shape[1], freq='h')
    columns = ['P','ET','GW', 'RO', 'theta_40cm', 'theta_60cm', 'theta_100cm']
    functions = 'mean', 'std', 'min', 'max', 'p5', 'p25', 'p75', 'p95'
    def do(s,data,axis):
         if s[0] == 'p':
             return np.percentile(data, int(s[1:]), axis=axis)
         else:
             return getattr(np, s)(data, axis=axis)
    df = pd.concat([
        pd.DataFrame(
            data=do(func, sim, axis=0),
            index=sim_index,
            columns=[c + '_' + func for c in columns]
        )
        for func in functions
    ], axis=1)
    return df

def run(code):
    with get_data(code) as data:
        threshold, _ = get_threshold(code, data)
        dotty_plot(code, data,threshold).savefig(f'fig/fig_{code}_dotty.png')
        plot_equifinality(code, data, threshold).savefig(f'fig/fig_{code}_balance.png')
        plot_validation(code, data, threshold).savefig(f'fig/fig_{code}_validation.png')
        plot_moisture(code, data, threshold).savefig(f'fig/fig_{code}_moisture.png')
        plot_boundary(code, data, threshold).savefig(f'fig/fig_{code}_boundary.png')
        result(code, data, threshold).to_excel(f'results_{code}.xlsx')  
        return get_tabdata(code, data, threshold)


def run_parallel(codes):
    from multiprocessing import Pool
    table = []
    i = 0
    with Pool(50) as pool:
        for code, tabdata in pool.imap_unordered(run, codes):
            i+=1
            print(f'{i:4}/{len(codes)} : {code}')
            table += tabdata
    pd.DataFrame(table).to_excel('run2.xlsx')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write('Usage: gbh_plot.py [11111|codes.txt]')
    elif os.path.exists(sys.argv[1]):
        run_parallel([l.strip() for l in open(sys.argv[1]) if '#' not in l])
    else:
        code = sys.argv[1]
        if '#' not in code:
           run(code)

    
