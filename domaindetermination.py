#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2023.06.20 V8 Voting & Variable Threshold Model - Patrick Kelley

def determinedomains(tempura,t,x,y,codefolderpath,df_com,LipidNames,mdnum,saveveryfig, randperc, raftcholdpsmpercthres, nonraftphospercthres,number_of_lipids):
    #################################### PARAMETERS ##########################################################
    # PLOT PARAMETERS
    # size of circle points
    pointsize = 150
    # border linewidth of circle points
    pointlnwth = 3

    # Same as Smoothing
    avgnumlipidinwindow = 11
    fracwindowslide = 3 # must be odd number

    # Phospholipid name
    phospholipid = ''
    for lip in LipidNames:
        if lip[-2:] == 'PC' or lip[-2:] == 'PE':
            if phospholipid == '':
                phospholipid = lip
            else:
                phospholipid = phospholipid+'-'+lip

    #########################################################################################################
    #################################### Import functions ###################################################
    #########################################################################################################
    import os, time, math, csv, copy
    import pandas as pd
    import numpy as np
    from functions import fillrin, fill2rin, threebythreebox, find_nearest, beads_read_gro, timelog, leafletdomaindetermine, leafletdomainareas
    # plots
    from functions import leafletscatterplot, leafletdomainplot
    import itertools
    from scipy.interpolate import griddata

    import sys
    sys.setrecursionlimit(6000)
    from sys import argv, stdout

    #########################################################################################################
    # FOLDER to store all the images and data
    picpath = codefolderpath+'/DomainDetermination/'
    isExist = os.path.exists(picpath)
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(picpath)
        stdout.write("New analysis directory is created! \n")

    # CSV datafile
    datapath = codefolderpath+'/DomainDetermination/'+'DomainDetermination.csv'  
    # Check if csv datafile exists, for phospholipid in phospholipidlist:if so, what was the last time entry & resume, otherwise start from scratch
    isExist = os.path.exists(datapath)
    if isExist:
        alldata = pd.read_csv(datapath)
        data = alldata['time[ps]'].values.tolist()
        # FIND NEAREST TIME IN SIMULATION TO TIMESNAP
        try:
            timesnap = int(data[-1])
            tempidx, temp = find_nearest(tempura,timesnap)
            saveveryfigindx,saveveryfigtemp = find_nearest(np.array(saveveryfig),timesnap) 
            if tempura[tempidx] != tempura[-1]:
                tempura = tempura[tempidx+1:]
                saveveryfig = saveveryfig[saveveryfigindx-1:]
            else:
                stdout.write('\nDomain Determination Data Directory is Complete. All domains have all been determined. Continuing...\n')
                return
        except (IndexError, ValueError):
            stdout.write('Domain Determination Data File is Empty. Deleting what is there.\n\n')
            os.remove(datapath)

    # Dataframe of raft/intermediate/nonraft domains, lipid ids inside each domain, and domain area and count
    df_columns = ['time[ps]', 'leaflet', 'comptime[s]', 
            'boxsize', 'windowsize', 'slidelength',
            'raftids', 'interids', 'nonraftids',
            'raftdomains', 'interdomains', 'nonraftdomains',
            'raftarea', 'interarea', 'nonraftarea', 
            'raftcnt', 'intercnt', 'nonraftcnt', 'totalarea']
    df_raft = pd.DataFrame(columns=df_columns)



    # subplt_n = 0 ### JOURNAL PLOT
    ##########################################################################################################
    ##########################################################################################################
    ################################## MAIN PROGRAM ##########################################################
    ##########################################################################################################
    ##########################################################################################################
    for temp in tempura:
        # IMAGES: to save or not to save:
        if temp in saveveryfig:
            figsave = True 
        else:
            figsave = False    

        totalcomptime_start = time.time()

        timeindex = t.index(temp)
        x_box = x[timeindex] #nm 
        y_box = y[timeindex] #nm 

        # Divide number_of_lipids by 2 because half in upper leaflet, half in lower leaflet
        xywind = math.sqrt((x_box*y_box)/((number_of_lipids)/2)*avgnumlipidinwindow)
        xyslide = xywind/fracwindowslide
        xwindow,ywindow = xywind,xywind
        xstep,ystep = xyslide,xyslide

        # if static grid (sliding full length of window); else sliding fraction of window
        if xstep == xwindow and ystep == ywindow:
            xsteplist = np.arange(0,x_box,xstep)
            ysteplist = np.arange(0,y_box,ystep)
        else:
            xsteplist = np.arange(0,x_box+xwindow/2,xstep)
            ysteplist = np.arange(0,y_box+ywindow/2,ystep)

        ##########################################################################################################
        ################################ UPPER LEAFLET ###########################################################
        ##########################################################################################################
        
        leaf = 'upper'
        uppcomptime_start = time.time()
        

        uppgridbox_x,uppgridbox_y,upperpoints,uppxi,uppyi,uppboxintplot,uppboxint,upperxydframe, uppraftids,uppintermedids,uppnonraftids = leafletdomaindetermine(df_com,leaf,temp, xstep,ystep,xwindow,ywindow,xsteplist,ysteplist, x_box,y_box, raftcholdpsmpercthres, nonraftphospercthres)

        # UPPER LEAFLET COMPUTATIONAL TIME
        uppcomptime_end = time.time()
        uppcomptime = uppcomptime_end - uppcomptime_start

        ##########################################################################################################
        ################################ LOWER LEAFLET ###########################################################
        ##########################################################################################################

        leaf = 'lower'
        lowcomptime_start = time.time()

        lowgridbox_x,lowgridbox_y,lowerpoints,lowxi,lowyi,lowboxintplot,lowboxint,lowerxydframe, lowraftids,lowintermedids,lownonraftids = leafletdomaindetermine(df_com,leaf,temp, xstep,ystep,xwindow,ywindow,xsteplist,ysteplist, x_box,y_box, raftcholdpsmpercthres, nonraftphospercthres)
        
        # LOWER LEAFLET COMPUTATIONAL TIME
        lowcomptime_end = time.time()
        lowcomptime = lowcomptime_end - lowcomptime_start

        ##########################################################################################################
        ##########################################################################################################
        ###################################### LEAFLET PLOTS #####################################################
        ##########################################################################################################
        ##########################################################################################################

        ##########################################################################################################
        #################################### LEAFLET COLOR CODE ##################################################
        ##########################################################################################################
        phosphocolor = ['darkblue', 'violet', 'navy']
        phosphos = phospholipid.split('-')
        ###### UPPER LEAFLET PLOT COLOR OF PHOSPHOLIPID, CHOLESTEROL, and SPHINGOMYELIN ##############
        upperspeed = ['orange']*len(upperxydframe['lipid'])
        for l,lip in enumerate(upperxydframe['lipid']):
            # DPSM
            if lip == 'DPSM':
                upperspeed[l] = 'red'
            # CHOL
            elif lip == 'CHOL':
                upperspeed[l] = 'white'
            # ATOC
            elif lip == 'ATOC':
                upperspeed[l] = '#FFBF00'
            # PHOSPHOLIPID
            elif lip == 'PUPC':
                upperspeed[l] = 'blue'
            else:
                for p, phospho in enumerate(phosphos):
                    if lip == phospho:
                        upperspeed[l] = phosphocolor[p]

        ###### LOWER LEAFLET PLOT COLOR OF PHOSPHOLIPID, CHOLESTEROL, and SPHINGOMYELIN ##############
        lowerspeed = ['orange']*len(lowerxydframe['lipid'])
        for l, lip in enumerate(lowerxydframe['lipid']):
            # DPSM
            if lip == 'DPSM':
                lowerspeed[l] = 'red'
            # CHOL
            elif lip == 'CHOL':
                lowerspeed[l] = 'white'
            # ATOC
# Colors: Orange - 'orange', Amber - '#FFBF00', Yellow - 'yellow', Pink - 'pink', Green - 'green', Black - 'black'
            elif lip == 'ATOC':
                lowerspeed[l] = '#FFBF00'
            # PHOSPHOLIPID
            elif lip == 'PUPC':
                lowerspeed[l] = 'blue'
            else:
                for p, phospho in enumerate(phosphos):
                    if lip == phospho:
                        lowerspeed[l] = phosphocolor[p]
        ##########################################################################################################
        # UPPER/LOWER LEAFLET RAFT/NONRAFT/INTERMEDIATE MAP
        ##########################################################################################################
        leafletscatterplot('upper',picpath,phospholipid,LipidNames,temp,randperc,figsave,upperpoints,upperspeed,pointlnwth,pointsize,x_box,y_box,xwindow,ywindow,xstep,ystep)
        leafletscatterplot('lower',picpath,phospholipid,LipidNames,temp,randperc,figsave,lowerpoints,lowerspeed,pointlnwth,pointsize,x_box,y_box,xwindow,ywindow,xstep,ystep)
        leafletdomainplot('upper',picpath,phospholipid,LipidNames,temp,randperc,figsave,uppxi,uppyi,uppboxintplot,upperpoints,upperspeed,pointlnwth,pointsize,x_box,y_box,xwindow,ywindow,xstep,ystep)
        leafletdomainplot('lower',picpath,phospholipid,LipidNames,temp,randperc,figsave,lowxi,lowyi,lowboxintplot,lowerpoints,lowerspeed,pointlnwth,pointsize,x_box,y_box,xwindow,ywindow,xstep,ystep)
        ##########################################################################################################
        # JOURNAL LOWER LEAFLET RAFT/NONRAFT/INTERMEDIATE MAP
        ##########################################################################################################
    #     subplt_n += 1    
    #     plt.style.use('default')
    #     plt.figure(figsize=(5*4, 10*4))
    #     plt.rcParams["figure.autolayout"] = True
    #     plt.rcParams.update({'font.size': 22})
    #     leafdomainplot(len(tempura),subplt_n,picpath,phospholipid,temp,figsave,lowxi,lowyi,lowboxintplot,lowerpoints,lowerspeed,pointlnwth,pointsize,randperc,xwindow,xstep,ywindow,ystep)

        ##########################################################################################################
        ############################# UPPER/LOWER LEAFLET DOMAIN AREAS ###########################################
        ##########################################################################################################
        # UPPER DOMAINS' AREAS
        ##########################################################################################
        uppraftdomain,uppinterdomain,uppnonraftdomain,uppraftA,uppinterA,uppnonraftA,uppraftAcount,uppinterAcount,uppnonraftAcount = leafletdomainareas(uppgridbox_x, uppgridbox_y, uppboxint)
        ##########################################################################################
        # LOWER DOMAINS' AREAS
        ##########################################################################################
        lowraftdomain,lowinterdomain,lownonraftdomain,lowraftA,lowinterA,lownonraftA,lowraftAcount,lowinterAcount,lownonraftAcount = leafletdomainareas(lowgridbox_x, lowgridbox_y, lowboxint)

        
        ##########################################################################################################
        ##########################################################################################################
        ####################################### SAVE DATA ########################################################
        ##########################################################################################################
        ##########################################################################################################

        ############################ SAVE UPPER RAFT IDS and AREA ################################################
        new_row = {'time[ps]':temp, 'leaflet':'upper', 'comptime[s]':uppcomptime, 
                    'boxsize':x_box, 'windowsize':xwindow, 'slidelength':xstep, 
                    'raftids':uppraftids, 'interids':uppintermedids, 'nonraftids':uppnonraftids,
                    'raftdomains': uppraftdomain, 'interdomains': uppinterdomain, 'nonraftdomains': uppnonraftdomain, 
                    'raftarea':uppraftA, 'interarea':uppinterA, 'nonraftarea':uppnonraftA, 
                    'raftcnt':uppraftAcount, 'intercnt':uppinterAcount, 'nonraftcnt':uppnonraftAcount, 
                    'totalarea':(uppraftA+uppinterA+uppnonraftA)}
        new_df = pd.DataFrame([new_row])
        df_raft = pd.concat([df_raft, new_df], axis=0, ignore_index=True)
        new_df = new_df[df_columns]
        # Check if csv datafile exists, if so, what was the last time entry & resume, otherwise start from scratch
        isExist = os.path.exists(datapath)
        if isExist:
            # append data frame to CSV file
            new_df.to_csv(datapath, mode='a', index=False, header=False)
        elif not isExist:
            # create the new data file
            with open(datapath, 'w') as csvfile:
                # using csv.writer method from CSV package
                write = csv.writer(csvfile)
                write.writerow(df_columns)
                # Close the file object
                csvfile.close()
            # append data frame to new data file
            new_df.to_csv(datapath, mode='a', index=False, header=False)

        ############################ SAVE LOWER RAFT IDS and AREA ################################################
        new_row = {'time[ps]':temp, 'leaflet':'lower', 'comptime[s]':lowcomptime,
                    'boxsize':x_box, 'windowsize':xwindow, 'slidelength':xstep, 
                    'raftids':lowraftids, 'interids':lowintermedids, 'nonraftids':lownonraftids,
                    'raftdomains': lowraftdomain, 'interdomains': lowinterdomain, 'nonraftdomains': lownonraftdomain,
                    'raftarea':lowraftA, 'interarea':lowinterA, 'nonraftarea':lownonraftA, 
                    'raftcnt':lowraftAcount, 'intercnt':lowinterAcount, 'nonraftcnt':lownonraftAcount, 
                    'totalarea':(lowraftA+lowinterA+lownonraftA)}
        new_df = pd.DataFrame([new_row])
        new_df = new_df[df_columns]
        df_raft = pd.concat([df_raft, new_df], axis=0, ignore_index=True)
        # Check if csv datafile exists, if so, what was the last time entry & resume, otherwise start from scratch
        isExist = os.path.exists(datapath)
        if isExist:
            # append data frame to CSV file
            new_df.to_csv(datapath, mode='a', index=False, header=False)
        elif not isExist:
            # create the new data file
            with open(datapath, 'w') as csvfile:
                # using csv.writer method from CSV package
                write = csv.writer(csvfile)
                write.writerow(df_columns)
                # Close the file object
                csvfile.close()
            # append data frame to new data file
            new_df.to_csv(datapath, mode='a', index=False, header=False)
        ##########################################################################################################
        totalcomptime_end = time.time()
        totalcomptime = totalcomptime_end - totalcomptime_start
        stdout.write('\n%s, time: %i, comptime: %0.2f'%(phospholipid,temp,totalcomptime))
        #################################### SAVE TO TIMELOG #####################################################
        timelog('dd',codefolderpath,phospholipid,mdnum,temp, totalcomptime_end, totalcomptime_start)
