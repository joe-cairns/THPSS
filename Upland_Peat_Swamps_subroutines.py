# python 3.8

# this code contains the subroutines used for the WRL Upland Peat Swamp daily water balance model
# do not make any changes to this code.
# this code needs to be in the same folder as the run code to work

import math
from math import tanh
import numpy as np
import pandas as pd
import datetime

def s_curves1(t, x4):
    """
        Unit hydrograph ordinates for UH1 derived from S-curves.
    """

    if t <= 0:
        return 0
    elif t < x4:
        return (t/x4)**2.5
    else: # t >= x4
        return 1


def s_curves2(t, x4):
    """
        Unit hydrograph ordinates for UH2 derived from S-curves.
    """

    if t <= 0:
        return 0
    elif t < x4:
        return 0.5*(t/x4)**2.5
    elif t < 2*x4:
        return 1 - 0.5*(2 - t/x4)**2.5
    else: # t >= x4
        return 1

def gr4j(precip, potential_evap, params, states = None, return_state = False):
    """
    #nbase code from https://github.com/amacd31/gr4j/blob/master/gr4j/gr4j.py.  has been adapted
        Generated simulated streamflow for given rainfall and potential evaporation.
        :param precip: Catchment average rainfall.
        :type precip: array(float)
        :param potential_evap: Catchment average potential evapotranspiration.
        :type potential_evap: array(float)
        :param params: X parameters for the model.
        :type params: dictionary with keys X1, X2, X3, X4
        :param states: Optional initial state values.
        :type states: Dictionary with optional keys 'production_store', 'routing_store'.
        :param return_state: If true returns a dictionary containing 'production_store' and 'routing_store'. Default: False.
        :type return_state: boolean
        :return: Array of simulated streamflow.
    """

    X1 = params['X1']
    X2 = params['X2']
    X3 = params['X3']
    X4 = params['X4']
    
    if states is None:
        states = {}

    nUH1 = int(math.ceil(X4))
    nUH2 = int(math.ceil(2.0*X4))

    uh1_ordinates = [0] * nUH1
    uh2_ordinates = [0] * nUH2

    UH1 = states.get('UH1', [0] * nUH1)
    UH2 = states.get('UH2', [0] * nUH2)

    for t in range(1, nUH1 + 1):
        uh1_ordinates[t - 1] = s_curves1(t, X4) - s_curves1(t-1, X4)

    for t in range(1, nUH2 + 1):
        uh2_ordinates[t - 1] = s_curves2(t, X4) - s_curves2(t-1, X4)

    production_store = states.get('production_store', 0) # S
    routing_store = states.get('routing_store', 0) # R
    
    qsim = []
    xsim=[]

    for P, E in zip(precip, potential_evap):

        if P > E:
            net_evap = 0
            scaled_net_precip = (P - E)/X1
            
            if scaled_net_precip > 13:
                scaled_net_precip = 13.
            tanh_scaled_net_precip = tanh(scaled_net_precip)
            reservoir_production = (X1 * (1 - (production_store/X1)**2) * tanh_scaled_net_precip) / (1 + production_store/X1 * tanh_scaled_net_precip)

            routing_pattern = P-E-reservoir_production
        else:
            
            scaled_net_evap = (E - P)/X1
            
            if scaled_net_evap > 13:
                scaled_net_evap = 13.
            tanh_scaled_net_evap = tanh(scaled_net_evap)
            
            ps_div_x1 = (2 - production_store/X1) * tanh_scaled_net_evap
            net_evap = production_store * (ps_div_x1) / (1 + (1 - production_store/X1) * tanh_scaled_net_evap)

            reservoir_production = 0
            routing_pattern = 0

        production_store = production_store - net_evap + reservoir_production
        
        groundwater_exchange = X2 * (routing_store / X3)**3.5

        #added to try matching Margot's model
        production_store = production_store
        percolation = production_store* (1 - 1 / (1 + (production_store/X1*4/9)**4)**0.25)
        
        routing_pattern = routing_pattern + (percolation)
        production_store = production_store - percolation
        
       
        for i in range(0, len(UH1) - 1):
            UH1[i] = UH1[i+1] + uh1_ordinates[i]*routing_pattern*0.9
        UH1[-1] = uh1_ordinates[-1] * routing_pattern*0.9

        for j in range(0, len(UH2) - 1):
            UH2[j] = UH2[j+1] + uh2_ordinates[j]*routing_pattern*0.1
        UH2[-1] = uh2_ordinates[-1] * routing_pattern*0.1
       
        routing_store = max(0, routing_store + UH1[0]+groundwater_exchange)
        

        R2 = routing_store / (1 + (routing_store / X3)**4)**0.25
        QR = routing_store - R2
        routing_store = R2
        QD = max(0, UH2[0]+groundwater_exchange)
        Q = QR + QD

        qsim.append(Q)#-groundwater_exchange)
        #qsim.append(-groundwater_exchange)
        xsim.append(groundwater_exchange)
    if return_state:
        return qsim, {
            'production_store': production_store,
            'routing_store': routing_store,
            'UH1': UH1,
            'UH2': UH2,
        }
    else:
        return qsim

       
def model2(params,gr4j_params_dict,stores,rain,ET,dates=[],allout=True,twoout=False): #allout=True for full output, false for just flow, twoout = True for flow and storage

    # run GR4J    
    flow, out_states=gr4j(rain, ET*gr4j_params_dict["Kco2"], gr4j_params_dict, states=stores, return_state=True)
    
    # Qin --> flow in from catchment (m3/day)
    # convert flow from mm/day to m3/day, by multiplying by the catchment area, less the swamp area
    Qin=np.array(flow)*(params["Ac"]-params["As"])/1000
    
    # dir_rain --> direct rain over the swamp area (m3/day)
    dir_rain=rain*params["As"]/1000
    
    # Qint is catchment rain plus the direct rainfall
    Qint = Qin + dir_rain
    
    # GWloss = groundwater loss 
    # in m3/day
    #GWloss = [params["Gwloss"]*params["As"] for i in rain]
    
    # ET_adj --> adjusted Et for crop coefficient, multiplied by swamp area to convert to m3/day
    ET_adj = ET*params["Kco"]*params["As"]/1000
    
    # calculate storages 
    surface_storage = params["SS"]
    effective_storage = params["SE"]
    total_storage = params["ST"]
    dead_storage = total_storage - surface_storage - effective_storage 
    
    #intial condition only water in the dead storage region, no percolation, ET in timestep 1
    St = total_storage 
    Se = effective_storage
    Sd = dead_storage
    perc1 = 0
    perc2 = 0
    Eloss = 0
    GW = 0
    Qout=[]
    
    eff_stor=[]
    tot_stor = []
    dead_stor = []
    sur_stor = []
    percK =[]
    percN =[]
    hor=[]
    sat=[]
    GWloss=[]
    AET=[]
    for Q, E, R in zip(Qint, ET_adj, rain): 
        if St<dead_storage:
            E=E*(St/dead_storage)#**2
            #E=0
        # calculate new total storage (St)
        if St - perc1 - perc2 + GW + Q - E  > 0:
            Q = Q 
            St = St - perc1 - perc2 + GW + Q - E

        else:
            St = 0
            
                
        #depending on the amount in storage, fill the different buckets --> Sd --> dead stoarge, Sa --> active storage, St --> total stoage, and runoff (excess rain goes straight to runoff)
        if St <= dead_storage:
            Sd = St
            Sa = 0
            Ss = 0
            sat_runoff = 0
        elif St < dead_storage + effective_storage:
            Sd = dead_storage
            Sa = St - Sd
            Ss = 0
            sat_runoff = 0
        elif St < total_storage:
            Sd = dead_storage
            Sa = effective_storage
            Ss = St - Sd - Sa
            sat_runoff = 0
        else:
            Sd = dead_storage
            Sa = effective_storage
            Ss = surface_storage
            sat_runoff = St - total_storage
            St = total_storage

        #percolation is proporation to current swamp storage 
        perc1 = min(params["Ke"]*params["Slope"]*(Sa/effective_storage)*Sa,Sa)
        perc2 = min(params["Ks"]*params["Slope"]*(Ss/surface_storage)*Ss,Ss)

        GW = -max(params["GW"]*(St/total_storage),0.6)*(params["As"]/1000)
        
        if twoout==True:
            tot_stor.append(St)
        if allout==True:
            eff_stor.append(Sa)
            tot_stor.append(St)
            dead_stor.append(Sd)
            sur_stor.append(Ss)
            percK.append(perc1)
            percN.append(perc2)
            sat.append(sat_runoff)
            GWloss.append(-GW)
            AET.append(E)#*St/total_storage)

        # percolation + any surface water runoff from excess rain equals total flow        
        Qout.append( sat_runoff + perc1 + perc2)   
    
    if allout == True:
        Qout=pd.DataFrame.from_dict({"Date":dates, "AET":AET, "GR4J out":Qin, "Inputs":Qint, "GW loss":GWloss, "Dead Storage":dead_stor, "Active Storage":eff_stor, "Surface Storage":sur_stor, "Total Storage":tot_stor, "Surface perc":percN, "Active storage perc": percK, "Saturation runoff":sat, "Modelled Flow":Qout,})
        return Qout
    elif twoout == True:
        return [Qout, tot_stor]
    else: return Qout

#objective functions        
def sqrtNSE(yo, ym):
    yo = np.array(yo)
    ym = np.array(ym)   
    #calculate Nash -Sutcliffe of Efficiency
    mask = np.array(~np.isnan(yo) & ~np.isnan(ym))
    y1 = yo[mask] 
    y2 = ym[mask]  
    y1mean = np.mean(np.power(y1,0.5))
    nse = 1 - np.sum(np.power(np.subtract(np.power(y2,0.5),np.power(y1,0.5)),2)) / np.sum(np.power(np.subtract(np.power(y1,0.5),np.power(y1mean,0.5)),2))
    return nse

def NSE(yo, ym):
    #calculate Nash -Sutcliffe of Efficiency
    yo = np.array(yo)
    ym = np.array(ym)
    yo = yo[:len(ym)]
    mask = np.array(~np.isnan(yo) & ~np.isnan(ym))
    y1 = yo[mask] 
    y2 = ym[mask]  
    y1mean = np.mean(y1)
    nse = 1 - np.sum(np.power(np.subtract(y2,y1),2)) / np.sum(np.power(np.subtract(y1,y1mean),2))
    return nse

def threeD(Q):
    Q3d = np.zeros(shape=len(Q)//3)
    for x in range(len(Q)//3):
        i = x*3
        if Q[i] != np.nan and Q[i+1] != np.nan and Q[i+2] != np.nan:
            Q3d[x]=(Q[i]+Q[i+1]+Q[i+2])
        else: Q3d[x]=np.nan
    return Q3d

def kge(eva, sim, return_all=False):
    """
    Kling-Gupta Efficiency
    Corresponding paper: 
    Gupta, Kling, Yilmaz, Martinez, 2009, Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling
    output:
        kge: Kling-Gupta Efficiency
    optional_output:
        cc: correlation 
        alpha: ratio of the standard deviation
        beta: ratio of the mean
    """
    eva = np.array(eva)
    sim = np.array(sim)
    mask = np.array(~np.isnan(eva) & ~np.isnan(sim))
    evaluation = eva[mask] 
    simulation = sim[mask] 
    if len(evaluation) == len(simulation):
        cc = np.corrcoef(evaluation, simulation)[0, 1]
        alpha = np.std(simulation) / np.std(evaluation)
        beta = np.sum(simulation) / np.sum(evaluation)
        kge = 1 - np.sqrt((cc - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
        if return_all:
            return kge, cc, alpha, beta
        else:
            return kge
    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan

def sqrtkge(eva, sim, return_all=False):
    """
    Kling-Gupta Efficiency
    Corresponding paper: 
    Gupta, Kling, Yilmaz, Martinez, 2009, Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling
    output:
        kge: Kling-Gupta Efficiency
    optional_output:
        cc: correlation 
        alpha: ratio of the standard deviation
        beta: ratio of the mean
    """
    eva = np.array(eva)
    sim = np.array(sim)
    mask = np.array(~np.isnan(eva) & ~np.isnan(sim))
    evaluation = np.power(eva[mask],0.5)
    simulation = np.power(sim[mask],0.5)
    if len(evaluation) == len(simulation):
        cc = np.corrcoef(evaluation, simulation)[0, 1]
        alpha = np.std(simulation) / np.std(evaluation)
        beta = np.sum(simulation) / np.sum(evaluation)
        kge = 1 - np.sqrt((cc - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
        if return_all:
            return kge, cc, alpha, beta
        else:
            return kge
    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan

def logkge(eva, sim, return_all=False):
    """
    Kling-Gupta Efficiency
    Corresponding paper: 
    Gupta, Kling, Yilmaz, Martinez, 2009, Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling
    output:
        kge: Kling-Gupta Efficiency
    optional_output:
        cc: correlation 
        alpha: ratio of the standard deviation
        beta: ratio of the mean
    """
    eva = np.array(eva)
    sim = np.array(sim)
    mask = np.array(~np.isnan(eva) & ~np.isnan(sim))
    evaluation = np.log(eva[mask])
    simulation = np.log(sim[mask])
    if len(evaluation) == len(simulation):
        cc = np.corrcoef(evaluation, simulation)[0, 1]
        alpha = np.std(simulation) / np.std(evaluation)
        beta = np.sum(simulation) / np.sum(evaluation)
        kge = 1 - np.sqrt((cc - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
        if return_all:
            return kge, cc, alpha, beta
        else:
            return kge
    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan
    
def bias(Qo,Qm):
    Qobst = 0
    Qmodt = 0
    Qo,Qm =np.array(Qo),np.array(Qm)
    mask = np.array(~np.isnan(Qo) & ~np.isnan(Qm))
    Qo = Qo[mask] 
    Qm = Qm[mask] 
    for Qob, Qmo in zip(Qo,Qm):
        if np.isnan(Qob)==False and Qmo <= 516:
            Qobst += Qob
            Qmodt += Qmo
    try:
        bias = 100*(Qmodt-Qobst)/Qobst
    except: bias = 1000
    return bias

def mixNSE(O,M):
    Qo = O[0]
    So = O[1]
    Qm = M[0]
    Sm = M[1]
    return min(sqrtNSE(Qo,Qm),NSE(So,Sm))

def mixNSE3d(O,M):
    Qo = O[0]
    So = O[1]
    Qm = M[0]
    Sm = M[1]
    return min(sqrtNSE(threeD(Qo),threeD(Qm)),NSE(So,Sm))

def mixNSE2(O,M,stor_min=0.5,breakdown=False):
    Qo = O[0]
    So = O[1]
    Qm = M[0]
    Sm = M[1]
    Q516 = O[2]
    if breakdown == False:
        #return sqrtNSE(Qo,Qm)*(0.5)+NSE(threeD(Q516),threeD(Qm))*(0.5) - abs(bias(Q516,Qm))/20 - max(0,(0.65- NSE(So,Sm)))
        return sqrtNSE(Qo,Qm)*(0.5)+NSE(threeD(Q516),threeD(Qm))*(0.5) - abs(bias(Q516,Qm))/10 - max(0,(0.7- NSE(So,Sm)))
    else:
        return sqrtNSE(Qo,Qm)*(0.5), NSE(threeD(Q516),threeD(Qm))*0.5, -abs(bias(Q516,Qm))/40, -max(0,(0.65- NSE(So,Sm)))
