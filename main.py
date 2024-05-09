#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2023.06.20 V8 Voting & Variable Threshold Model - Patrick Kelley
# 2022.02.23 V7 generalize for all lipids & combine all analysis - Patrick Kelley
# 2021.06.12 V6 Domain Determination - Patrick Kelley
# 2021.28.07 V1 - Patrick
# 2021.31.05 raft-leaflet-order-paramater.py UPDATE V2 - Patrick Kelley
# 2021.11.08 CG CoM cos2 NN - Patrick Kelley
# UPDATE NEEDED: read box size from Structure folder gro files, end -> BOX SIZE function

##########################################################################################################
# !!!!!!!!!!!!!!!!!!! place path for all the .xtc, ~.edr~ and .tpr/.gro files !!!!!!!!! !!!!!!!!!!!!!!!!!!
##########################################################################################################
#datafolderpath = './run1.1/'
datafolderpath = './run1/'

##########################################################################################################
##########################################################################################################
###################### IMPORT PYTHON PACKAGES ############################################################
##########################################################################################################
##########################################################################################################

import math
import os
import time
import numpy as np
import pandas as pd
from sys import argv, stdout

##########################################################################################################
#################### Define Parameters of Simulation #####################################################
##########################################################################################################
# gmx check -f #file#.xtc -> Timestep (ps) = 500 = 0.02*25000
dt = 0.02
nstout= 25000

# bilayer normal
orientation_of_bilayer_normal = [0, 0, 1]
# number of lipids
number_of_lipids = 1010
# lipid list
LipidNames = ['ATOC','PUPC','DPSM','CHOL']
# Backbone Beads CoM instead of Global CoM
BackBoneBeadCondition = True

# Time jumps ~10ns
desiredvmdstep = 10000
vmdstep = round(desiredvmdstep/(nstout*dt))
# Save a figure every desirefigsavesteps ~ 100ns
figsaveinterval = 100000 # every 100ns

# gmx check #frames*Timestep (ps) = 6000300 = 13,334 * 450ps
prodruntime = 5000000 #5micros

codefolderpath = './'
boxpath = codefolderpath+'boxsize.csv'

############################## DEFINE VARIABLE THRESHOLD VALUES ##########################################
randperc = 0.25000 # 25% Minimum Probability (i.e. 1 in 4 chance of seeing this randomly)
maxwindowlipid = 50

# PATRICK'S POPC/ATOC SIMULATION Leaflet 
lipidchar = ['A','B','C','D'] # ATOC, PUPC,DPSM,CHOL
# Raft-like lipids: [CHOL/DPSM]-rich ['C','D']
raftlike = ['C','D']
# Nonraft-like lipids: PL[PDPC]-rich ['B']
nonraftlike = ['B']

# Multiply by two to prevent N frXom being less than m in binomial distribution in choose r lipids > available lipid types
# otherwise math.factorial(m)*math.factorial(N-m) would be a negative factorial !ERROR!
lipidnumber = [50*2, 320*2, 320*2, 320*2] # ATOC, PDPC, DPSM, CHOL

from functions import variablethresholdvalues
raftcholdpsmpercthres, nonraftphospercthres = variablethresholdvalues(maxwindowlipid, lipidchar,lipidnumber,randperc,raftlike, nonraftlike, codefolderpath)
##########################################################################################################

# INPUT FILE LIST
xtc_list = sorted([f for f in os.listdir(datafolderpath) if os.path.isfile(os.path.join(datafolderpath, f)) and f.endswith('.xtc')])
tpr_list = sorted([f for f in os.listdir(datafolderpath) if os.path.isfile(os.path.join(datafolderpath, f)) and f.endswith('.tpr')])
if not tpr_list:
    tpr_list = sorted([f for f in os.listdir(datafolderpath) if os.path.isfile(os.path.join(datafolderpath, f)) and f.endswith('.gro')])
mdnumlist = range(len(xtc_list))


#########################################################################################################

# (Normalized) Orientation of Bilayer Normal
binorm = math.sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
for i in range(3):
    orientation_of_bilayer_normal[i] /= binorm
traj_skip = vmdstep

df_com = pd.DataFrame(data={"time": [], "lipid": [], "id": [], "resid_start": [], "resid_end": [], "leaflet": [],  "x": [], "y": [], "z": []})
t,x,y,z = [],[],[],[]

##########################################################################################################
########################################## MAIN COMPUTATION ##############################################
##########################################################################################################
# IMPORT center of mass and [cos(theta)]^2 function
from centerofmass import com
# IMPORT domain determination function
from domaindetermination import determinedomains
# IMPORT figure saving at time increments 
from functions import savefigurestimesincrements

for mdnum in mdnumlist:
    start_frame = time.time()
    
    trajfile = datafolderpath + xtc_list[mdnum]
    tprfile = datafolderpath + tpr_list[mdnum]
    ############################## TIME IN FILE CORRECTION for jump steps  ####################################
    initial_time0 = (mdnum+0)*prodruntime
    final_time0 = (mdnum+1)*prodruntime
    initial_time = int(math.ceil(initial_time0/(dt*nstout))*(dt*nstout))
    final_time = range(int(math.ceil(initial_time0/(dt*nstout))*(dt*nstout)),final_time0,int((vmdstep)*(dt)*(nstout)))[-1]
    stdout.write("\n"+"*"*(40))
    stdout.write("\nStarting on MD: %i \t %i-%i ps\n"%(mdnum,initial_time,final_time))
    stdout.write("*"*(40)+"\n")
    ########################################################################################################## 
    # Create folder for all analysis
    pathname = codefolderpath+'Analysis_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps/'
    # Create folder for structure files
    structpathname = codefolderpath+'Structure_Files_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps/'
    
    ##########################################################################################################
    #################################  CoM, [cos(theta)]^2 & BoxSize #########################################
    ##########################################################################################################
    df_com_new = com(initial_time,vmdstep,dt,nstout,final_time,prodruntime, orientation_of_bilayer_normal,BackBoneBeadCondition,LipidNames,trajfile,tprfile,pathname,traj_skip,codefolderpath,mdnum)
    df_com = pd.concat([df_com,df_com_new], ignore_index=True)
    ##########################################################################################################

    ##########################################################################################################
    ###################################  Box Size  ###########################################################
    ##########################################################################################################
    df_boximport = pd.read_csv(boxpath)
    t = df_boximport['t[ps]'].loc[np.logical_or((df_boximport['t[ps]'] <= initial_time),(df_boximport['t[ps]'] <= final_time))].values.tolist()
    x = df_boximport['x'].loc[np.logical_or((df_boximport['t[ps]'] <= initial_time),(df_boximport['t[ps]'] <= final_time))].values.tolist()
    y = df_boximport['y'].loc[np.logical_or((df_boximport['t[ps]'] <= initial_time),(df_boximport['t[ps]'] <= final_time))].values.tolist()
    z = df_boximport['z'].loc[np.logical_or((df_boximport['t[ps]'] <= initial_time),(df_boximport['t[ps]'] <= final_time))].values.tolist()
    
    ##########################################################################################################        

    ##########################################################################################################
    ################################# Domain Determination ###################################################
    ##########################################################################################################
    tempura = np.unique(df_com['time']) # List of time[ps]
    # Save Figure Every Interval Amount
    saveveryfig = savefigurestimesincrements(tempura,figsaveinterval)
    # Voting & Variable Threshold Domain Determination
    determinedomains(tempura,t,x,y,codefolderpath,df_com,LipidNames,mdnum,saveveryfig, randperc, raftcholdpsmpercthres, nonraftphospercthres,number_of_lipids)

    ##########################################################################################################


##########################################################################################################
#################################### MAIN ANALYSIS #######################################################
##########################################################################################################   
##########################################################################################################
################################# Domain Order Parameters ################################################
##########################################################################################################
from domainorderparameter import determinedomainorderparameters
df_raft_op_dict, df_inter_op_dict, df_nonraft_op_dict = determinedomainorderparameters(codefolderpath, pathname, LipidNames, mdnumlist,dt,nstout,vmdstep,prodruntime)    
        
