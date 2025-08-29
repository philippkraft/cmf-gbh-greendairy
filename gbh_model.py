"""
"""
import cmf
import pandas as pd
import spotpy
import os
import numpy as np
import time

def load_weather(p: cmf.project, fn: str):
    """
    Load weather data from Excel file and transform it into meteo and
    rain stations

    """
    weather_df = pd.read_excel(fn, index_col=0)
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

def u(name: str, lower: float, upper:float, default:float=None,doc=None):
    "Creates a spotpy parameter with uniform distribution"
    return spotpy.parameter.Uniform(name, lower, upper, optguess=default, doc=doc)


class Structure:
    """
    Stellt die Modellstrukturtypen da
    """
    def __init__(self, modelcode):
        self.modelcode = str(modelcode)
        self.ET = int(self.modelcode[4])
        self.PERC = int(self.modelcode[3])
        self.CANOPY = bool(int(self.modelcode[2])-1)
        self.SURF = bool(int(self.modelcode[1])-1)
        self.MACRO = bool(int(self.modelcode[0])-1)

    def __str__(self):
        return self.modelcode
    



class Model:
    def __init__(self, structure_code: str):
        self.structure = Structure(structure_code)
        path = '.' # os.path.dirname(__file__)
        # self.data_fn = os.path.join(path, 'weatherstation.xlsx')
        self.soilmoisture = load_soilmoisture('moisture.xlsx')        
        cmf.set_parallel_threads(1)

        p = self.project = cmf.project()
        self.cell: cmf.Cell = p.NewCell(0,0,0, 1000)
        self.meteo, self.rain = load_weather(p, 'weather.xlsx')
        self.create_nodes()
        
        self.parameters = [
                u('soil_capacity',
                    10, 1000, default=500,
                    doc='Capacity of rooted zone in :math:`mm`'
                ),
                u('ETV1',
                    0.1, 0.9, default=0.7,
                    doc='Fraction of soil capacity where ET starts to be limited'
                ),
                u('ETV0', 0.0, 0.9, default=0.35, doc='Faction of ETV1 where wilting point is reached'),
            ]
        if self.structure.PERC == 3:
            self.parameters += [
                u('vg_alpha', 0.005, 0.1, default=0.01, doc='Van Genuchten alpha value'),
                u('vg_n', 1.1, 1.8, default=1.2, doc='Van Genuchten alpha value'),
                ]
        else:
            self.parameters += [
                u('Sfc', 0.1, 1.0, default=0.7, doc='Fraction of capacity at field capacity in mm/mm. field capacity = Sfc * Smax'),
                u('retention_time',
                    1, 150, default=10,
                    doc='Retention time in the soil water storage'
                )
            ]
        if self.structure.PERC == 3 or self.structure.SURF:
            self.parameters += [u('Ksat', 1, 200, default=50, doc='Ksat in mm/h')]
        
        if self.structure.MACRO:
            self.parameters += [u('infexp', 0.001, 2.0, default=0.1, doc='Infiltration exponent for BROOK90 macropore infiltration')]
        self.macropore_connections = []


        
        
    def create_nodes(self):
        self.groundwater : cmf.DirichletBoundary = self.project.NewOutlet('gw', 0,0,-2)
        self.runoff : cmf.DirichletBoundary = self.project.NewOutlet('runoff', 0,0,0)
        if self.structure.PERC == 1:
            self.cell.add_layer(1)
        elif self.structure.PERC == 2:
            for d in [0.1, 0.3, 0.5, 0.7, 1.0]:
                self.cell.add_layer(d)
        elif self.structure.PERC == 3:
            for d in np.arange(0.02, 1.01, 0.02):
                self.cell.add_layer(d, cmf.VanGenuchtenMualem())
        if self.structure.CANOPY:
            self.canopy_storage_evaporation()
        if self.structure.MACRO:
            self.macro = self.cell.add_storage('macropores')

            
    
    def canopy_storage_evaporation(self):  
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
 

    def set_soil_capacity(self, p) -> float:
        """
        Sets the upper soil capacity
        """
        for layer in self.cell.layers:
            # Parameterize soil water capacity
            if self.structure.PERC == 3:
                layer.soil.Phi = p.soil_capacity / self.cell.soildepth / 1000 
            else:
                layer.soil.porosity = p.soil_capacity / self.cell.soildepth / 1000 
            # Set initial value
            layer.wetness = 0.5
        return p.soil_capacity

        
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
    

    def parameterize_et(self, p):
        self.cell.set_uptakestress(cmf.ContentStress(p.ETV1, p.ETV1 * p.ETV0))
        if self.structure.ET == 1:
            cmf.PriestleyTaylorET.use_for_cell(self.cell)
        else:
            cmf.PenmanMonteithET.use_for_cell(self.cell)
            
    def parameterize_perc(self, p):
        if self.structure.SURF or self.structure.PERC == 3:
            for layer in self.cell.layers:
                layer.soil.Ksat = p.Ksat * 24 / 1000

        
        if self.structure.PERC <= 2:
            cmf.LinearStorageConnection(
                self.cell.layers[-1], 
                self.groundwater,
                p.retention_time, 
                residual = p.Sfc * self.cell.layers[-1].get_capacity()
            )
            for upper, lower in zip(self.cell.layers[:-1], self.cell.layers[1:]):
                cmf.LinearStorageConnection(
                    upper, 
                    lower,
                    p.retention_time, 
                    residual = p.Sfc * upper.get_capacity()
                )
                
        elif self.structure.PERC == 3:
            cmf.FreeDrainagePercolation(self.cell.layers[-1], self.groundwater)
            for upper, lower in zip(self.cell.layers[:-1], self.cell.layers[1:]):
                cmf.Richards(upper, lower)
            for layer in self.cell.layers:
                layer.soil.alpha = p.vg_alpha
                layer.soil.n = p.vg_n
                layer.soil.fit_w0(w1=1.01, Psi_p=1.0)

    def parameterize_surf(self, p):
        if self.structure.SURF:
            cmf.ConceptualInfiltration(
                self.cell.layers[0],
                self.cell.surfacewater, 
                0.95
                )
            cmf.WaterbalanceFlux(self.cell.surfacewater, self.runoff)
            
        if self.structure.MACRO:
            # Get the flux fraction that should go to the deeper layers (=1-fraction for layer 0)
            infiltration_to_deeper = 1 - self.inf_fraction(self.cell.layers[0], p.infexp)
            pfr0 = cmf.PartitionFluxRoute(self.cell.surfacewater, self.cell.layers[0], self.macro, infiltration_to_deeper)
            pfr = [
                cmf.PartitionFluxRoute(self.cell.layers[0], self.macro, l, self.inf_fraction(l, p.infexp)/pfr0.fraction)
                for l in self.cell.layers[1:]
                ]
            return [pfr0] + pfr



    def parameterize(self, p):
        self.set_soil_capacity(p)
        self.parameterize_et(p)
        self.parameterize_perc(p)
        self.macropore_connections = self.parameterize_surf(p)
            

            
        
    def AMBAV_rss(self, JT: int):
        #Simulation der Reifephase, Rss = Faktor 
        JTs = 90
        JTz = 120
        JTmax=180  #Jahrestag des maximalen Wassersverbrauch 
        JTe = 220
        
        rss = rssmin = 46 #Wert aus AMBAV Dokumentation -> ist auch rss min? 
        
        if JT>JTe or JT<JTs:
            lai=0.3

        elif JT>JTmax:
            lai = 3.8
            rss = rssmin + (
                60 * (
                    (JT-JTmax)/(JTmax+45-JTmax)
                    )
                ) + (              
                150 * (
                    (JT-JTmax)/(JTmax+45-JTmax)
                    )**3
                )
        elif JT<JTz:
            lai = 0.3 + (0.7 - 0.3)* (JT-JTs)/(JTz-JTs)
        elif JT<=JTmax:
            lai = 0.7 + (3.8-0.7) * (JT-JTz)/(JTmax-JTz)
        
        return lai, rss
        
   
    def Kc_WinterWeizen(self, JT: int):
        JTs = 90
        JTmax=180  #Jahrestag des maximalen Wassersverbrauch 
        JTe = 220
        if JT>JTe or JT<JTs:
            return 0.4
        elif JT<JTmax:
            return 0.4 + (1.15 - 0.4)* (JT-JTs)/(JTmax-JTs)
        else:
            return 0.4 + (1.15 - 0.4)* (JTe-JT)/(JTe-JTmax)


    def simulation(self, p, verbose=False):
        
        self.parameterize(p)
        solver = cmf.CVodeBanded(self.project, 1e-9)
        result = []
        tstart = time.time()
        for t in solver.run(self.meteo.T.begin, self.meteo.T.end, cmf.h):
            # Kill to slow running models
            if time.time() - tstart > 300:
                return np.NaN * np.zeros((len(self.meteo.T), 6))
            ET = self.cell.transpiration(t) + self.cell.evaporation(t)
            GW = self.groundwater(t)
            P = -self.cell.rain_source(t)
            RO = self.runoff(t)
            if self.structure.ET==2:
                self.cell.vegetation.Kc = self.Kc_WinterWeizen(t.DOY())
            if self.structure.ET==3:
                self.cell.vegetation.LAI, self.cell.vegetation.rss = self.AMBAV_rss(t.DOY())
                self.cell.vegetation.height = self.cell.vegetation.LAI / 3.8
            if verbose and not (solver.t % (cmf.day)).is_not_0():
                V = self.cell.vegetation
                print(f'{t:%Y-%m-%d %H:%M} : P={P:0.3f} ET={ET:0.3f} GW={GW:0.3f} RO={RO:0.3f}, LAI={V.LAI:0.2f}, Kc={V.Kc:0.2f}')
            if self.structure.PERC == 1:
                outputlayers = [self.cell.layers[0]] * 2
            elif self.structure.PERC == 2:
                outputlayers = [self.cell.layers[2], self.cell.layers[3]]
            elif self.structure.PERC == 3:
                outputlayers = [self.cell.layers[19], self.cell.layers[29]]
            result.append([P, ET, GW, RO] + [l.theta for l in outputlayers])
            
        return np.array(result)
    def simulation_to_df(self, simulation):
        sim_index = self.meteo.T.to_pandas().index
        return pd.DataFrame(
            data=simulation, 
            index=sim_index, 
            columns = ['P','ET','GW', 'RO', 'theta_40cm', 'theta_60cm']
        )
        
    
    def evaluation(self):
        
        calib = self.soilmoisture[0]
        return np.array(calib[['d40cm_mean', 'd60cm_mean']])

    def mean_on_perc1(self, data):
        if self.structure.PERC == 1:
            return np.array(data).mean(1)
        else:
            return np.array(data)
        
    def objectivefunction(self, simulation, evaluation, params=None):
        if simulation is None:
            return [np.NaN] * 6
        sim_df = self.simulation_to_df(simulation)
        def rmse(eval, sim):
            return np.sqrt(np.nanmean((eval-sim)**2))

        compare = [
                [
                    self.mean_on_perc1(sm[['d40cm_mean', 'd60cm_mean']]), 
                    self.mean_on_perc1(sim_df[sim_df.columns[-2:]].loc[sm.index])
                ]
                for sm in self.soilmoisture
        ]
        def objf(ev, sim):
            try:
                return rmse(ev, sim)
            except:
                return np.NaN

        return [
            objf(np.array(cp[0]), np.array(cp[1]))
            for cp in compare
        ] + list(simulation[-365*24:,:4].mean(0) * 365)

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('Usage: python gbh_model.py <model_code>')
    else:
        code = sys.argv[1]
        model = Model(code)
        p = spotpy.parameter.create_set(model, 'optguess')
        model.parameterize(p)
        print(cmf.describe(model.project))
        sim = model.simulation(p, verbose=True)
        sim_df = model.simulation_to_df(sim)
        objf = model.objectivefunction(sim, model.evaluation())
        print(f'RMSE40cm={objf[0]:0.3f}, RMSE60cm={objf[1]:0.3f}, '
              f'P={objf[2]:0.3f}, ET={objf[3]:0.3f}, GW={objf[4]:0.3f}, RO={objf[5]:0.3f}')
        print(sim_df.tail())


        # print(model.objectivefunction(model.simulation(spotpy.parameter.create_set(model, 'random')), model.evaluation()))