"""
"""
import cmf
import pandas as pd
import spotpy
import os
import numpy as np
import time
import typing
import re

def load_weather(p: cmf.project, fn: str):
    """
    Load weather data from Excel file and transform it into meteo and
    rain stations

    """
    weather_df = pd.read_excel(fn, index_col=0)
    weather_df.columns = [
        col[:col.index(' at ')].lower().replace(' ','_').replace('.','') 
        for col in weather_df.columns
    ]
    weather_df.loc[pd.isna(weather_df.global_radiation), 'global_radiation'] = 0.0
    weather_df.ffill(inplace=True)
    meteo: cmf.MeteoStation = p.meteo_stations.add_station(
        'gbh', (0,0,0), 50.39, 8.24,1.0, 
        pd.to_datetime(weather_df.index[0]), cmf.h)
    meteo.T = cmf.timeseries.from_pandas(weather_df.air_temperature)
    meteo.Tmax = meteo.T.floating_max(cmf.day)
    meteo.Tmin = meteo.T.floating_min(cmf.day)
    meteo.InstrumentHeight = 2
    meteo.Rs = cmf.timeseries.from_pandas(weather_df.global_radiation * 86400 * 1e-6)
    meteo.Sunshine = cmf.timeseries.from_pandas(weather_df.sunshine_fraction)
    meteo.rHmean = cmf.timeseries.from_pandas(weather_df.rel_humidity)
    meteo.Windspeed = cmf.timeseries.from_pandas(weather_df.wind_speed)
    
    rain = p.rainfall_stations.add(
        'gbh', 
        cmf.timeseries.from_pandas(weather_df.rainfall * 24),
        (0,0,0)
    )
    p.use_nearest_meteo()
    p.use_nearest_rainfall()

    return meteo, rain

def load_soilmoisture(fn: str):
    sm_df = pd.read_excel(fn, index_col=0)

    # Remove problematic data    
    sm_calib = pd.DataFrame(index=sm_df[sm_df.index<'2024-05-01'].index, columns=['d40cm_mean','d40cm_std','d60cm_mean','d60cm_std'])
    sm_valid = pd.DataFrame(index=sm_df[sm_df.index>'2024-05-01'].index, columns=['d40cm_mean','d40cm_std','d60cm_mean','d60cm_std'])
    
    sm_calib.d40cm_mean = sm_df[sm_df.columns[::2]].mean(1)
    sm_calib.d40cm_std = sm_df[sm_df.columns[::2]].std(1)
    sm_calib.d60cm_mean = sm_df[sm_df.columns[1::2]].mean(1)
    sm_calib.d60cm_std = sm_df[sm_df.columns[1::2]].std(1)
    sm_valid.d40cm_mean = sm_df[sm_df.columns[::2]].mean(1)
    sm_valid.d40cm_std = sm_df[sm_df.columns[::2]].std(1)
    sm_valid.d60cm_mean = sm_df[sm_df.columns[1::2]].mean(1)
    sm_valid.d60cm_std = sm_df[sm_df.columns[1::2]].std(1)
    
    return sm_calib, sm_valid

def load_kc(fn: str):
    kc_df = pd.read_excel(fn, sheet_name='Kc_Verlauf', index_col='date')
    kci = kc_df['Kc'].resample('d').interpolate()
    rdi = kc_df['RootDepth80'].resample('d').interpolate()
    return cmf.timeseries.from_pandas(kci), cmf.timeseries.from_pandas(rdi)

def u(name: str, lower: float, upper:float, default:typing.Optional[float]=None,doc=None):
    "Creates a spotpy parameter with uniform distribution"
    return spotpy.parameter.Uniform(name, lower, upper, optguess=default, doc=doc)

def ug(name: str, mean: float, sd:float, default:typing.Optional[float]=None, doc=None):
    "Creates a spotpy parameter with uniform distribution around mean +/- 3*sd"
    return spotpy.parameter.Uniform(name, mean - 3 * sd, mean + 3 * sd, optguess=default or mean, doc=doc)

class Parameters:
    DW = ug('DW', mean=18, sd=2, doc='Dead water in % of soil (> 4.2 pF)')
    aFC = ug('aFC', mean=14.5, sd=1.2, doc='Plant available field capacity % of soil (4.2 - 2.5 pF)')
    AC = ug('AC', mean=9.6, sd=2.5, doc='Air capacity in % of soil (< 2.5 pF)')
    Ksat = u('Ksat', lower=0.05, upper=1, doc='Saturated hydraulic conductivity in m/d')
    percexp = u('percexp', 1.0, 5.0, default=2.0, doc='Exponent for the power law percolation')
    infexp = u('infexp', 0, 0.5, default=0, doc='Infiltration exponent for BROOK90 macropore infiltration')
    
    @classmethod
    def as_list(cls):
        return [p for p in vars(cls).values() if isinstance(p, spotpy.parameter.Base)]
    



class Model:
    def __init__(self):
        self.soilmoisture = load_soilmoisture('moisture.xlsx') 
        self.Kc, self.RootDepth = load_kc('Kc.xlsx')      
        cmf.set_parallel_threads(1)

        p = self.project = cmf.project()
        self.cell: cmf.Cell = p.NewCell(0,0,0, 1000)
        self.meteo, self.rain = load_weather(p, 'weather-2025.xlsx')
        self.cell.vegetation.fraction_at_rootdepth = 0.8
        self.create_nodes()
        
        self.parameters = Parameters.as_list()
        self.macropore_connections = []

   
        
    def create_nodes(self):
        self.groundwater : cmf.DirichletBoundary = self.project.NewOutlet('gw', 0,0,-2)
        self.runoff : cmf.DirichletBoundary = self.project.NewOutlet('runoff', 0,0,0)
        for d in [0.1, 0.3, 0.5, 0.7, 0.9, 1.1]:
            self.cell.add_layer(d)
        # self.canopy_storage_evaporation()
        self.macro = self.cell.add_storage('macropores')

            
    
    def canopy_storage_evaporation(self):  
        raise ValueError("Canopy storage not used in GBH model")
        # Make c.canopy a water storage
        self.cell.add_storage('Canopy','C')
        # Split the rainfall from the rain source (RS) between 
        # intercepted rainfall (RS->canopy) and throughfall (RS-surface)
        cmf.Rainfall(self.cell.canopy, self.cell, False, True) # RS->canopy, only intercepted rain
        cmf.Rainfall(self.cell.surfacewater, self.cell, True, False) # RS->surface, only throughfall
        # Use an overflow mechanism, eg. the famous Rutter-Interception Model
        cmf.RutterInterception(self.cell.canopy, self.cell.surfacewater, self.cell) 
        # And now the evaporation from the wet canopy (using a classical Penman equation)
        cmf.CanopyStorageEvaporation(self.cell.canopy, self.cell.evaporation, self.cell)
 

    def set_soil_capacity(self, p: Parameters) -> float:
        """
        Sets the upper soil capacity
        """
        for layer in self.cell.layers:
            # Parameterize soil water capacity
            if layer.upper_boundary < 0.1:
                AC = p.AC + 4
            elif layer.upper_boundary < 0.3:
                AC = p.AC + 2.5
            else:
                AC = p.AC 
            porosity = (p.DW + p.aFC + AC) / 100 # in m3/m3
            layer.soil.porosity = porosity
            # Set initial value
            layer.wetness = 0.5

        
    def inf_fraction(self, layer: cmf.SoilLayer, infexp: float, max_depth=None):
        """
        Calculates the infiltration fraction for soil layer as in BROOK90
        
        see: http://www.ecoshift.net/brook/wat.html#INFPAR
        

        Parameters
        ----------
        layer : cmf.SoilLayer
            The cmf layer to calculate the infiltration fraction.
        infexp : float
            The infiltration exponent.
            
        max_depth : float, optional
            The maximum depth for infiltration (IDEPTH in BROOK90). If None,
            use the depth of the deepest soil layer in cmf
            The default is None.

        Returns
        -------
        float
            The infiltration fraction for layer

        """
        l = layer.lower_boundary
        u = layer.upper_boundary
        t = max_depth or layer.cell.soildepth
        if infexp==0:  # if infexp is 0, all water goes to the first layer
            return float(u==0)
        elif l > t:
            return 0.0
        return (l/t)**infexp - (u/t)**infexp    
    

    def parameterize_et(self, p: Parameters):
        ETV1 = (p.DW + p.aFC / 2) / 100
        self.cell.set_uptakestress(cmf.ContentStress(ETV1, p.DW / 100))
        cmf.PenmanMonteithET.use_for_cell(self.cell)

    def make_matrix_connection(self, upper: cmf.SoilLayer, lower: cmf.WaterStorage, p: Parameters):
        C = upper.get_capacity()
        Sfc = (p.aFC + p.DW) / 100 * upper.thickness * upper.cell.area  # in m³
        return cmf.PowerLawConnection(
            upper, 
            lower,
            Q0 = p.Ksat * upper.cell.area,  # in mm/d
            V0 = C - Sfc,
            beta = p.percexp,
            residual = Sfc
        )
            
    def parameterize_perc(self, p: Parameters):
        for layer in self.cell.layers:
            layer.soil.Ksat = p.Ksat  # in m/d
        self.perc_connections = [
            self.make_matrix_connection(upper, lower, p)
            for upper, lower in zip(self.cell.layers, list(self.cell.layers[1:]) + [self.groundwater])
        ]

    def parameterize_surf(self, p):
        cmf.ConceptualInfiltration(
            self.cell.layers[0],
            self.cell.surfacewater, 
            0.95
            )
        cmf.WaterbalanceFlux(self.cell.surfacewater, self.runoff)
            
        # Get the flux fraction that should go to the deeper layers (=1-fraction for layer 0)
        if hasattr(p, 'infexp') and p.infexp > 0:
            infiltration_to_deeper = 1 - self.inf_fraction(self.cell.layers[0], p.infexp)
            pfr0 = cmf.PartitionFluxRoute(self.cell.surfacewater, self.cell.layers[0], self.macro, infiltration_to_deeper)
            pfr = [
                cmf.PartitionFluxRoute(self.cell.layers[0], self.macro, l, self.inf_fraction(l, p.infexp)/pfr0.fraction)
                for l in self.cell.layers[1:]
                ]
            return [pfr0] + pfr
        else:
            return []


    def parameterize(self, p):
        self.set_soil_capacity(p)
        self.parameterize_et(p)
        self.parameterize_perc(p)
        self.macropore_connections = self.parameterize_surf(p)

    def result(self, t=None):
        """
        Returns the actual model state at time t as a numpy array
        """
        t = t or self.meteo.T.begin
        return np.array([
            -self.cell.rain_source(t),
            self.cell.transpiration(t) + self.cell.evaporation(t),
            self.groundwater(t),
            self.runoff(t),
        ] + [l.theta for l in self.cell.layers[-3:]])
    
    def simulation_shape(self):
        return len(self.meteo.T), len(self.result())
            
   
    def simulation(self, p, verbose=False):
        
        self.parameterize(p)
        solver = cmf.CVodeDense(self.project, 1e-9)
        result = []
        tstart = time.time()
        for t in solver.run(self.meteo.T.begin, self.meteo.T.end, cmf.h):
            # Kill too slow running models
            if time.time() - tstart > 300:
                return np.nan * np.zeros((len(self.meteo.T), 6))
            ET = self.cell.transpiration(t) + self.cell.evaporation(t)
            GW = self.groundwater(t)
            P = -self.cell.rain_source(t)
            RO = self.runoff(t)
            self.cell.vegetation.Kc = self.Kc[t]
            self.cell.vegetation.RootDepth = self.RootDepth[t]
            if verbose and not (solver.t % (cmf.day)).is_not_0():
                V = self.cell.vegetation
                print(f'{t:%Y-%m-%d %H:%M} : P={P:0.3f} ET={ET:0.3f} GW={GW:0.3f} RO={RO:0.3f}, LAI={V.LAI:0.2f}, Kc={V.Kc:0.2f}, RootD={V.RootDepth:0.2f}')
            result.append(self.result(t))
            
        return np.array(result)
    
    def simulation_to_df(self, simulation):
        sim_index = self.meteo.T.to_pandas().index
        return pd.DataFrame(
            data=simulation, 
            index=sim_index, 
            columns = ['P','ET','GW', 'RO', 'theta_40cm', 'theta_60cm', 'theta_100cm']
        )
        
    
    def evaluation(self):
        
        calib = self.soilmoisture[0]
        return np.array(calib[['d40cm_mean', 'd60cm_mean']])
        
    def objectivefunction(self, simulation, evaluation, params=None):
        if simulation is None:
            return [np.nan] * 6
        sim_df = self.simulation_to_df(simulation)
        def rmse(eval, sim):
            return np.sqrt(np.nanmean((eval-sim)**2))

        compare = [
                [
                    sm[['d40cm_mean', 'd60cm_mean']], 
                    sim_df[sim_df.columns[-2:]].loc[sm.index]
                ]
                for sm in self.soilmoisture
        ]
        def objf(ev, sim):
            try:
                return rmse(ev, sim)
            except:
                return np.nan

        return [
            objf(np.array(cp[0]), np.array(cp[1]))
            for cp in compare
        ] + list(simulation[-3*365*24:,:4].mean(0) * 365)

if __name__ == '__main__':
    import sys
    model = Model()
    p = spotpy.parameter.create_set(model, 'optguess')
    model.parameterize(p)
    print(cmf.describe(model.project))
    sim = model.simulation(p, verbose=True)
    sim_df = model.simulation_to_df(sim)
    objf = model.objectivefunction(sim, model.evaluation())
    print(f'RMSE40cm={objf[0]:0.3f}, RMSE60cm={objf[1]:0.3f}, '
            f'P={objf[2]:0.3f}, ET={objf[3]:0.3f}, GW={objf[4]:0.3f}, RO={objf[5]:0.3f}')
    print(sim_df.tail())
