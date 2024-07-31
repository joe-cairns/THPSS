# python 3.8 

# This is the code that is used to run the Upland Peat Swamp daily water balance model developed by WRL
# the Upland_Peat_Swamps_subroutines.py should be in the same folder to run


import pandas as pd
import Upland_Peat_Swamps_subroutines as THPSS_subroutines
import matplotlib.pyplot as plt, mpld3
import numpy as np
import datetime
import matplotlib
import os
from functools import reduce

#input file  # THIS IS THE ONLY THING THAT NEEDS TO BE CHANGED
# PUT INPUT FILE INTO THE SAME FOLDER AS THE CODE, AND REPLACE NAME HERE

data_file= ".xlsx"

#DO NOT CHANGE BELOW THIS LINE
##############################################################################################################################################
#model inputs
rain=pd.read_excel(data_file, sheet_name="Rain", parse_dates=["Date"])
ET = pd.read_excel(data_file, sheet_name="ET", parse_dates=["Date"])
Qobs = pd.read_excel(data_file, sheet_name="Qobs", parse_dates=["Date"])
storage=pd.read_excel(data_file, sheet_name="Storage", parse_dates=["Date"])

gp= pd.read_excel(data_file, sheet_name="gr4j_params", index_col="Parameter")
mp= pd.read_excel(data_file, sheet_name="swamp_params", index_col="Parameter")
    
#paramater dictionaries
params = {"Location":mp["Value"]["Location name"],"SS":mp["Value"]["Ss"], "SE":mp["Value"]["Se"],"ST":mp["Value"]["St"], "Ke":mp["Value"]["Ke"], "Ks":mp["Value"]["Ks"], "As":mp["Value"]["As"], "Ac":mp["Value"]["Ac"], "Slope":mp["Value"]["Slope"], "Kco":mp["Value"]["Kco"], "In":mp["Value"]["In"],"GW":mp["Value"]['GW']}
gr4j_params_dict={"X1":gp["Value"]["X1"], "X2":gp["Value"]["X2"], "X3":gp["Value"]["X3"], "X4":gp["Value"]["X4"], "Kco2":gp["Value"]["Kco"]}
stores={"production_store":gp["Value"]["X1"]*0.3, "routing_store":gp["Value"]["X3"]*0.3}

# make the directories for the outputs
out_file = os.path.join("Outputs", params["Location"], '%s output.xlsx' % params["Location"])
if not os.path.exists(os.path.dirname(out_file)):
    os.makedirs(os.path.dirname(out_file))

out_graphs = os.path.join("Outputs", params["Location"], '%s output.png' % params["Location"])
if not os.path.exists(os.path.dirname(out_graphs)):
    os.makedirs(os.path.dirname(out_graphs))

# use the rain and the ET to get the dates as daily rain and ET are needed
# dates only include times when both are available
rain_date=rain["Date"]
ET_date=ET["Date"]

start=max(min(rain_date), min(ET_date))
end=min(max(rain_date), max(ET_date))

rain=rain[(rain["Date"]>=str(start.date()))&(rain["Date"]<=str(end.date()))]
ET=ET[(ET["Date"]>=str(start.date()))&(ET["Date"]<=str(end.date()))]

# now check for flow and storage in the same date range
storage=storage[(storage["Date"]>=str(start.date()))&(storage["Date"]<=str(end.date()))]
Qobs=Qobs[(Qobs["Date"]>=str(start.date()))&(Qobs["Date"]<=str(end.date()))]

# combine all the data into one
data_combined = reduce(lambda  left,right: pd.merge(left,right,on=['Date'], how='outer'), [rain, ET, storage, Qobs])
data_combined=data_combined[(data_combined["Date"]>=str(start.date()))&(data_combined["Date"]<=str(end.date()))]

#full model outputs
Qout = THPSS_subroutines.model2(params,gr4j_params_dict,stores,rain["Daily Rainfall"],ET["Daily ET"],allout=True,dates=rain["Date"])

# this adds secondary flows which aren't routed through the swmap
flow, states= THPSS_subroutines.gr4j(rain["Daily Rainfall"], ET["Daily ET"]*gr4j_params_dict["Kco2"], gr4j_params_dict, states=stores, return_state=True)

Qout["Secondary Catchment"]=np.array(flow)*mp["Value"]["A2"]/1000
Qout["Modelled Flow"]= Qout["Modelled Flow"] + Qout["Secondary Catchment"]                             

# collated data for the output spreadsheet
data_sheet = pd.DataFrame(data={
    "Date": data_combined['Date'],
    "Observed flow (m3/day)": data_combined['Observed flow (m3/day)'],
    "Modelled flow (m3/day)": Qout["Modelled Flow"],
    "Observed Storage (m3)": data_combined['Observed storage'],
    "Modelled Storage (m3)": Qout['Total Storage'],
    "AET": Qout["AET"],  # Include the actual column name used in Qout
    "GR4J out": Qout["GR4J out"],  # Include the actual column name used in Qout
    "Inputs": Qout["Inputs"],  # Include the actual column name used in Qout
    "GW loss": Qout["GW loss"],  # Include the actual column name used in Qout
    "Inactive Storage": Qout["Dead Storage"],  # Include the actual column name used in Qout
    "Active Storage": Qout["Active Storage"],  # Include the actual column name used in Qout
    "Surface Storage": Qout["Surface Storage"],  # Include the actual column name used in Qout
    "Surface perc": Qout["Surface perc"],  # Include the actual column name used in Qout
    "Active storage perc": Qout["Active storage perc"],  # Include the actual column name used in Qout
    "Saturation runoff": Qout["Saturation runoff"],  # Include the actual column name used in Qout
})


# Do some stats
if len(Qobs)<1:
    Q_provided='No'
    NSE_sqrt_flow=np.nan
    NSE_low_flow=np.nan
    NSE_flow=np.nan
else:
    Q_provided='Yes'
    
    #observed flows less than 1000
    Qobs1000 = data_sheet.copy()
    Qobs1000.loc[Qobs1000['Observed flow (m3/day)']>1000, "Observed flow (m3/day)"] = np.nan
    #flows for different obj. function calcs
    Qo = data_sheet["Observed flow (m3/day)"]
    Qm= data_sheet["Modelled flow (m3/day)"]


    NSE_sqrt_flow=THPSS_subroutines.sqrtNSE(Qo,Qm)
    NSE_low_flow=THPSS_subroutines.sqrtNSE(Qobs1000['Observed flow (m3/day)'],Qm)
    NSE_flow=THPSS_subroutines.sqrtNSE(Qo,Qm)

if len(storage)<1:
    #storage for different obj. function calcs
    S_provided = 'No'
    NSE_storage=np.nan
else:
    S_provided = 'Yes'
    Sm = data_sheet["Modelled Storage (m3)"]
    So = data_sheet["Observed Storage (m3)"]
    
    NSE_storage=THPSS_subroutines.sqrtNSE(So, Sm)

stats=pd.DataFrame({"Value":[Q_provided, S_provided, NSE_low_flow,NSE_flow,NSE_sqrt_flow, NSE_storage]}, index=["Observed flow data provided?","Observed storage data provided?", "NSE flow (<1,000m3/day)", "NSE flow","NSE square root of flow","NSE Storage"])    

# write the outputs to an excel file
with pd.ExcelWriter(out_file) as writer:  
    data_sheet.to_excel(writer, sheet_name="Model outputs", index=False)
    stats.to_excel(writer, sheet_name='Statistics')
    mp.to_excel(writer, sheet_name='Input Swamp_params')
    gp.to_excel(writer, sheet_name='Input GR4J_params')

#plots
font = { 'size'   : 8}

matplotlib.rc('font', **font)
fig, ax = plt.subplots(2, 1, figsize=(8,12), sharex=True)

if Q_provided=='Yes':
    data_sheet.plot(x="Date", y="Observed flow (m3/day)", ax=ax[0], label="Observed Discharge", color="darkorange")
data_sheet.plot(x="Date", y="Modelled flow (m3/day)", ax=ax[0], label="Modelled Discharge (intact)",linestyle="--", color="royalblue")

if S_provided=='Yes':
    data_sheet.plot(x="Date", y="Observed Storage (m3)", ax=ax[1], label="Observed Storage", color="seagreen")
data_sheet.plot(x="Date", y="Modelled Storage (m3)", ax=ax[1], label="Modelled Storage (intact)",linestyle="--", color="firebrick")
ax[0].xaxis.set_tick_params(which='both', labelbottom=True)
plt.savefig(out_graphs, dpi=300)