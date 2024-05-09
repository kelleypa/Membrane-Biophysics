#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2022.02.23 V7 generalize for all lipids & combine all analysis - Patrick Kelley
##########################################################################################################
######################### IMPORT/CALCULATE&SAVE ORDER PARAMETER  #########################################
##########################################################################################################

def determinedomainorderparameters(codefolderpath, pathname, LipidNames, mdnumlist,dt,nstout,vmdstep,prodruntime):
    import os
    import csv
    # Import writer class from csv module
    # from csv import writer
    import time
    import pandas as pd
    import math
    import numpy as np
    import copy
    from functions import find_csv_filenames
    
    import sys
    sys.setrecursionlimit(6000)
    from sys import argv, stdout
    
    from lipidinventory import lipidinventory
    df_beadmass,bond_list,beads_list,backbonebeads_list = lipidinventory(LipidNames)
    
    # Phospholipid name
    phospholipid = ''
    for lip in LipidNames:
        if lip[-2:] == 'PC' or lip[-2:] == 'PE':
            if phospholipid == '':
                phospholipid = lip
            else:
                phospholipid = phospholipid+'-'+lip
    leaflet = ["upper","lower"]      
    # LIPID IDS IN DOMAINS CSV datafile
    raftnonraftpathname = codefolderpath+'DomainDetermination/DomainDetermination.csv'  
    # Read RAFT, NONRAFT, AREAS, TOTALAREA CSV if it exists:
    if os.path.exists(raftnonraftpathname) == True:
        df_domain = pd.read_csv(raftnonraftpathname)

        # INDEX REMINDER (DPSM <-> PHOSPHOLIPID: VMD)

        # create dictionary of raft and nonraft order parameter (op)
        df_raft_op_dict,df_nonraft_op_dict, df_inter_op_dict = {},{},{}
        df_raft_op = pd.DataFrame(data=np.unique(df_domain['time[ps]']), columns=['time'])
        df_inter_op = pd.DataFrame(data=np.unique(df_domain['time[ps]']), columns=['time'])
        df_nonraft_op = pd.DataFrame(data=np.unique(df_domain['time[ps]']), columns=['time'])
        for leaff in leaflet:
            df_raft_op_dict[leaff] = {}
            df_inter_op_dict[leaff] = {}
            df_nonraft_op_dict[leaff] = {}
            for molecool in LipidNames:
                df_raft_op_dict[leaff][molecool] = df_raft_op.copy()
                df_inter_op_dict[leaff][molecool] = df_inter_op.copy()
                df_nonraft_op_dict[leaff][molecool]= df_nonraft_op.copy()

        #########################################################################################################

        # SUM RAFT/NONRAFT/INTERMEDIATE OR IMPORT

        for leaf in leaflet: 
            for mol, molecul in enumerate(LipidNames):
                opraftpathname = codefolderpath+'/DomainDetermination/'+phospholipid+ '_' + leaf + '_' + molecul +'_RaftOrderParameters.csv'  
                opinterpathname = codefolderpath+'/DomainDetermination/'+phospholipid+ '_' + leaf + '_' + molecul +'_MixedOrderParameters.csv'  
                opnonraftpathname = codefolderpath+'/DomainDetermination/'+phospholipid+ '_' + leaf + '_' + molecul +'_NonraftOrderParameters.csv'  
                # CHECK if OP CSV datafiles exists and if so, import
                opraftpathexists,opinterpathexists,opnonraftpathexists = os.path.exists(opraftpathname),os.path.exists(opinterpathname),os.path.exists(opnonraftpathname)
                if opraftpathexists == True and opinterpathexists == True and opnonraftpathexists == True:
                    df_raft_op_dict[leaf][molecul] = pd.read_csv(opraftpathname)  
                    #stdout.write("\n3/2[cos(theta)]^2-1/2 for %s %s - raft - exists. Importing. \n" % (leaf, molecul))

                    df_inter_op_dict[leaf][molecul] = pd.read_csv(opinterpathname)  
                    #stdout.write("\n3/2[cos(theta)]^2-1/2 for %s %s - mixed - exists. Importing. \n" % (leaf, molecul))

                    df_nonraft_op_dict[leaf][molecul] = pd.read_csv(opnonraftpathname)   
                    #stdout.write("\n3/2[cos(theta)]^2-1/2 for %s %s - nonraft - exists. Importing. \n" % (leaf, molecul))

                for mdnum in mdnumlist:
                    ############################## TIME IN FILE CORRECTION for jump steps  ####################################
                    initial_time0 = mdnum*prodruntime
                    final_time0 = (mdnum+1)*prodruntime
                    initial_time = int(math.ceil(initial_time0/(dt*nstout))*(dt*nstout))
                    final_time = range(int(math.ceil(initial_time0/(dt*nstout))*(dt*nstout)),final_time0,int((vmdstep)*(dt)*(nstout)))[-1]
                    ############################################################################################################
                    # if the 3/2[cos(theta)]^2-1/2 file exists, check if last time is in raft file
                    if opraftpathexists == True and opinterpathexists == True and opnonraftpathexists == True:
                        if final_time in df_raft_op_dict[leaf][molecul]['time'].values:
                            stdout.write("\n[cos(theta)]^2 for %s %s - MDNUM=%i has already been summed; moving to next production run.\n" % (leaf, molecul, mdnum))
                            continue
                    # else compute the 3/2[cos(theta)]^2-1/2 for the production run and append to file
                    else:
                        # Locate folder for theta analysis at specific production run [ie mdrun]
                        pathname = codefolderpath+'Analysis_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps/'

                        filenames = find_csv_filenames(pathname)

                        bonds = bond_list[mol].split()

                        for bond_name in bonds:
                            stdout.write("\n\tSumming 3/2[cos(theta)]^2-1/2 for %s %s - %s - MDNUM=%i.\n" % (leaf, molecul, bond_name, mdnum))
                            file1 = open(pathname+molecul+'_theta_'+bond_name+'.dat','r')
                            Lines = file1.readlines()
                            # separate spaces between string data points
                            strLines = [Lines[i].split(',') for i in range(0,len(Lines))]
                            # remove '\n' from lines
                            for i, sline in enumerate(strLines):
                                for j, svalue in enumerate(sline):
                                    strLines[i][j] = svalue.strip('\n')
                            # convert into arrays
                            strArray = np.array(strLines)
                            # convert empty places (ie lipid not in the leaflet) to np.nan
                            strArray[strArray == ''] = np.nan
                            # convert strings to floats
                            X = strArray[1:].astype(np.float64)
                            # create dataframe from array with columns defined by title line
                            df = pd.DataFrame(X, columns = strArray[0])
                            # all ids for the molecule
                            df_id = df.columns

                            for indtim, timmie in enumerate(np.unique(df['time'])):
                                if int(timmie) in np.unique(df_domain['time[ps]']):
                                    df_time = df.loc[df['time'] == timmie]
                                    df_raft_time = df_domain['raftids'].loc[np.logical_and(df_domain['time[ps]'] == timmie, df_domain['leaflet'] == leaf)].tolist()[0]
                                    # remove '][' and filter out ''
                                    df_raft_time = list(filter(None,map(str.strip, df_raft_time.strip('][').replace('"', '').split(','))))
                                    df_raft_time = [int(float(i)) for i in df_raft_time]
                                    df_inter_time = df_domain['interids'].loc[np.logical_and(df_domain['time[ps]'] == timmie, df_domain['leaflet'] == leaf)].tolist()[0]
                                    df_inter_time = list(map(str.strip, df_inter_time.strip('][').replace('"', '').split(',')))
                                    df_inter_time = [int(float(i)) for i in df_inter_time]
                                    df_nonraft_time = df_domain['nonraftids'].loc[np.logical_and(df_domain['time[ps]'] == timmie, df_domain['leaflet'] == leaf)].tolist()[0]
                                    # remove '][' and filter out ''
                                    df_nonraft_time = list(filter(None,map(str.strip, df_nonraft_time.strip('][').replace('"', '').split(','))))
                                    df_nonraft_time = [int(float(i)) for i in df_nonraft_time]
                                    raftiesum,raft_cnt = 0,0
                                    intersum,inter_cnt = 0,0
                                    nonraftiesum,nonraft_cnt = 0,0
                                    for ndx in df_id:
                                        # RAFT ORDER PARAMETER
                                        if ndx == 'time':
                                            continue
                                        elif int(ndx) in df_raft_time:
                                            raftiesum = raftiesum + 3/2*np.cos(df[str(ndx)].loc[df['time']==timmie].values[0])**2-1/2
                                            raft_cnt += 1
                                        # INTERMEDIATE ORDER PARAMETER
                                        if ndx == 'time':
                                            continue
                                        elif int(ndx) in df_inter_time:
                                            intersum = intersum + 3/2*np.cos(df[str(ndx)].loc[df['time']==timmie].values[0])**2-1/2
                                            inter_cnt += 1
                                        # NONRAFT ORDER PARAMETER
                                        if ndx == 'time':
                                            continue
                                        elif int(ndx) in df_nonraft_time:
                                            nonraftiesum = nonraftiesum + 3/2*np.cos(df[str(ndx)].loc[df['time']==timmie].values[0])**2-1/2
                                            nonraft_cnt += 1
                                    df_raft_op_dict[leaf][molecul].loc[df_raft_op_dict[leaf][molecul]['time']==timmie, 'N'] = raft_cnt
                                    df_raft_op_dict[leaf][molecul].loc[df_raft_op_dict[leaf][molecul]['time']==timmie, bond_name] = raftiesum
                                    df_inter_op_dict[leaf][molecul].loc[df_inter_op_dict[leaf][molecul]['time']==timmie, 'N'] = inter_cnt
                                    df_inter_op_dict[leaf][molecul].loc[df_inter_op_dict[leaf][molecul]['time']==timmie, bond_name] = intersum
                                    df_nonraft_op_dict[leaf][molecul].loc[df_nonraft_op_dict[leaf][molecul]['time']==timmie, 'N'] = nonraft_cnt
                                    df_nonraft_op_dict[leaf][molecul].loc[df_nonraft_op_dict[leaf][molecul]['time']==timmie, bond_name] = nonraftiesum
                        # SAVE TO CSV
                        # CSV datafile
                        for datapath in [opraftpathname, opinterpathname, opnonraftpathname]:
                            with open(datapath, mode='a', newline='') as csvfile:
                                if datapath == opraftpathname:
                                    df_raft_op_dict[leaf][molecul].to_csv(csvfile, header=(csvfile.tell()==0))
                                elif datapath == opinterpathname:
                                    df_inter_op_dict[leaf][molecul].to_csv(csvfile, header=(csvfile.tell()==0))
                                elif datapath == opnonraftpathname:
                                    df_nonraft_op_dict[leaf][molecul].to_csv(csvfile, header=(csvfile.tell()==0))
                            stdout.write('Domain determination order parameter data for %s leaflet, %s lipid of MDNUM=%i has been saved in %s.\n'%(leaf, molecul,mdnum,datapath))
                   
        stdout.write("\n\n Finished, please come again. %s\n" % (" "*56))
        return df_raft_op_dict, df_inter_op_dict, df_nonraft_op_dict
    else:
        stdout.write("Domain Determination Data missing. Terminating. %s\n" % (" "*56))
