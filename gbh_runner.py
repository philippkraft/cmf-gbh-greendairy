# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:19:51 2024

@author: philipp
"""

import pandas as pd
import spotpy
import numpy as np
import time
import datetime
import tables
import logging
from gbh_model import Model

logger = logging.getLogger(__name__)

class Runner:
    def get_table_def(self, params, sim_shape):
        """
        Returns a dict of column definitions using multidimensional
        hdf5 columns. Columns for parameters and likelihoods are atomic and resemble
        the csv datawriter. If ``save_sim=True``, the simulation array is saved as an array value in
        a single multidimensional table cell

        cf.: https://www.pytables.org/usersguide/tutorials.html#multidimensional-table-cells-and-automatic-sanity-checks

        """
        # Start position of parameter columns
        param_pos = 6
        # Start position of simulation columns
        sim_pos = param_pos + len(params)
        dtype=np.dtype('float32')
        names = ['RMSEcalib', 'RMSEvalid', 'P', 'ET', 'GW', 'RO'] + list(params.name)
        columns = {
            name: tables.Col.from_dtype(dtype, pos=i)
            for i, name in enumerate(names)
        }

        # Get the appropriate dtype for the n-d cell
        # (tables.Col.from_dtype does not take a shape parameter)
        sim_dtype = np.dtype((dtype, sim_shape))
        columns['simulation'] = tables.Col.from_dtype(
            sim_dtype, pos=sim_pos
        )
        logger.info(f"Created table definition with columns: {list(columns.keys())}")
        return columns

    def __init__(self, code):
        """
        Creates the runnner for a specific code, not much is happening here
        """
        self.code = code
        self.db = None
        self.table = None
        self.dbname = f'gbh_{code}'
        
        self.header = []

    def __enter__(self):
        """
        Enters a context block (with) and opens the hdf5-table for writing

        returns the Runner with an opened HDF5-Table
        """
        self.db = tables.open_file(f'{self.dbname}.h5', 'w', self.dbname)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Closes the table on exit of the with block
        """
        self.db.close()


    def sample(self, runs):
        """
        Samples n runs using the LHS method, rewritten from spotpy. Returns a generator generating a code, parameterset tuple.
        If the table is not yet created, it does so now-

        Specific to the gbh_model, eg. hardcoded shape of simulation result
        """
        model = Model()
        parameters = spotpy.parameter.get_parameters_array(model)
        parameterset = spotpy.parameter.ParameterSet(parameters)
        if not self.table:

            sim_shape = model.simulation_shape()
            columns = self.get_table_def(parameterset, sim_shape)
            self.table = self.db.create_table('/', 'data', description=columns)
            self.header = list(columns.keys())

        segment = 1 / runs
        matrix = np.empty((runs, len(parameters)))
        segments = np.arange(0, 1, segment)
        param_range = parameters['maxbound'] - parameters['minbound']
        for i in range(len(parameters)):
            matrix[:, i] = parameters['minbound'][i] + (segments + np.random.uniform(0, segment, runs)) * param_range[i]
            np.random.shuffle(matrix[:,i])
        
        for row in matrix:
            yield parameterset(*row)

    
    def save(self, objf, params, result):
        """
        Saves a single run to the HDF5 table
        """
        new_row = self.table.row
        for header, value in zip(self.header, list(objf) + list(params) + [result]):
            new_row[header] = value
        new_row.append()
        

def run_lhs(name=None, runs=None):
    """
    Runs the gbh_model[code] the given number of runs or as default
    - 10000 runs for single storage models
    - 2000 runs for tipping bucket models
    - 500 runs for Richard's equation models
    """
    from multiprocessing import Pool
    runs = runs or 2000
    name = name or str(runs)
    i = saved = 0
    min_obj = 1000
    tstart = time.time()
    start = datetime.datetime.now()
    with Pool(50) as pool:
        with Runner(name) as runner:
            for objective, params, result in pool.imap_unordered(run_tuple, runner.sample(runs)):
                try:
                    i += 1
                    if objective[0] < 0.05:
                        runner.save(objective, params, result)
                        saved += 1
                    min_obj = min(min_obj, objective[0])
                    print(f'{i:5d}/{runs:5d}', f'RMSEmin: {min_obj:7.4f}, RMSEnow: {objective[0]:7.4f}, {time.time()-tstart:6.0f}sec, {saved} runs saved')
                except KeyboardInterrupt:
                    break
    runtime = (time.time() - tstart) / 60
   

def run_tuple(parameterset):
    model = Model()
    try:
        result = model.simulation(parameterset)
        objf = model.objectivefunction(result, ...)
    except:
        result = np.zeros(model.simulation_shape()) * np.nan
        objf = [np.nan] * 6
    return objf, list(parameterset), result

def run(parameterset=None):
    model = Model()
    parameterset = parameterset or spotpy.parameter.create_set(model, 'random')
    tstart = time.time()
    result = model.simulation(parameterset, verbose=False)
    duration =time.time() - tstart
    return duration, model.objectivefunction(result, ...), parameterset



if __name__ == '__main__':
    import sys
    today = datetime.datetime.now().strftime('%Y-%m-%d')
    run_lhs(f'KA6-{today}', 10000)
    from gbh_plot import run as plot_run
    plot_run(f'KA6-{today}')

