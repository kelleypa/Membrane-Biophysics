#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2022.02.23 V7 generalize for all lipids & combine all analysis - Patrick Kelley


##########################################################################################################
##########################################################################################################
############################## calculate CoM and theta  #########################################
##########################################################################################################
##########################################################################################################   
def com(initial_time,vmdstep,dt,nstout,final_time,prodruntime, orientation_of_bilayer_normal,BackBoneBeadCondition,LipidNames,trajfile,tprfile,pathname,traj_skip,codefolderpath,mdnum):
    from functions import beads_read_gro, boxsize_read_gro
    import os
    from os import remove
    import time
    import pandas as pd
    from sys import argv, stdout
    import subprocess
    import math
    from functions import timelog, saveordervalues
    ##########################################################################################################
    comstrFile = pathname+'CoM_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps.dat'
    isExist = os.path.exists(comstrFile)
    if BackBoneBeadCondition is True:
        comstrFile = pathname+'BBB_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps.dat'
        isExist = os.path.exists(comstrFile)
    
    # Check whether the specified (analysis) pathname exists or not
    isanalysisExist = os.path.exists(pathname)
    if not isanalysisExist:
        # Create a new directory because it does not exist 
        os.makedirs(pathname)
        stdout.write("New analysis directory is created! \n")
    ##########################################################################################################

    if not isExist:
        start_total = time.time()
        
        ##########################################################################################################
        stdout.write("Starting center of mass calculation & theta calculation \n\n")
        
        ##########################################################################################################
        ################## Output all frame using trjconv  (aka timesnaps) #######################################
        ##########################################################################################################
        # Create folder for structure files
        structpathname = './Structure_Files_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps/'
        # Check whether the specified pathname exists or not
        isStructExist = os.path.exists(structpathname)
        if not isStructExist:
            # Create a new directory because it does not exist 
            os.makedirs(structpathname)
            stdout.write("New structure (.gro) directory is created! \n\n")
        # Output all frame using trjconv 
        dir = os.listdir(structpathname)
        if len(dir) == 0:
            stdout.write("Output all coordinate files \n")
            #command = "gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %s/frame_dump_.gro > /dev/null" % (trajfile, tprfile, initial_time, final_time, traj_skip,structpathname)
            try:
                command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %sframe_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time, final_time, traj_skip, structpathname)
                #command = "echo %s | gmx trjconv -f %s -s %s  -o %s.whole.xtc -pbc whole > /dev/null" % ("System", trajfile, tprfile, trajfile)
                #command = "echo %s | gmx trjconv -f %s.whole.xtc -s %s -b %i -e %i -sep -skip %i -pbc whole -o %sframe_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time, final_time, traj_skip, structpathname)
                print(command)
                subprocess.check_output(command, shell=True)
                #subprocess.call(command, shell=True)
                continuationflag = True
                initial_time_restart = initial_time
                final_time_restart = final_time
            except subprocess.CalledProcessError as e:
                ############################## TIME IN FILE CORRECTION for time starting over at 0 microseconds ####################################
                initial_time00 = 0*prodruntime
                final_time00 = (0+1)*prodruntime
                initial_time_restart = int(math.ceil(initial_time00/(dt*nstout))*(dt*nstout))
                final_time_restart = range(int(math.ceil(initial_time00/(dt*nstout))*(dt*nstout)),final_time00,int((vmdstep)*(dt)*(nstout)))[-1]
                command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %sframe_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time_restart, final_time_restart, traj_skip,structpathname)
                print(command)
                subprocess.call(command, shell=True)
                continuationflag = False

        else:
            stdout.write("\nStructure folder is not empty. Should we run gromacs 'gmx trjconv' to find .gro structure files? Yes:anykey , No:n \n")
            #choice = input("> ")  
            choice = "n"
            if choice == "n":
                stdout.write("Continuing... \n\n")
                continuationflag = False
            else: 
                stdout.write("Output all coordinate files \n")
                try:
                    command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %sframe_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time, final_time, traj_skip,structpathname)
                    print(command)
                    subprocess.check_output(command, shell=True)
                    subprocess.call(command, shell=True)
                    continuationflag = True
                    initial_time_restart = initial_time
                    final_time_restart = final_time
                except subprocess.CalledProcessError as e:
                    ############################## TIME IN FILE CORRECTION for time starting over at 0 microseconds ####################################
                    initial_time00 = 0*prodruntime
                    final_time00 = (0+1)*prodruntime
                    initial_time_restart = int(math.ceil(initial_time00/(dt*nstout))*(dt*nstout))
                    final_time_restart = range(int(math.ceil(initial_time00/(dt*nstout))*(dt*nstout)),final_time00,int((vmdstep)*(dt)*(nstout)))[-1]
                    command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %sframe_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time_restart, final_time_restart, traj_skip,structpathname)
                    print(command)
                    subprocess.call(command, shell=True)
                    continuationflag = False
        
        ##########################################################################################################
        ################################### TOP-BOTTOM LEAFLETS  #################################################
        ##########################################################################################################
        leaflet = ["upper","lower"]     
        from lipidinventory import lipidinventory
        df_beadmass,bond_list,beads_list,backbonebeads_list = lipidinventory(LipidNames)
        
        ##########################################################################################################
        ##########################################################################################################
        ################## For each dumped frame (aka timesnap) ##################################################
        ##########################################################################################################
        file_count = 0
        com_count = 0
        while True:
            ###################### CLEAR/DEFINE THETA (bond), DX, DY, DZ PARAMETERs ########################################################
            # Nested Dictionary of Bonds dictionaries inside Lipids dictionaries
            df_bond,df_dx,df_dy,df_dz = {}, {},{},{}
            ###################### INITIALIZE CENTER OF MASS PARAMETERs #########################################################
            
            df_com = pd.DataFrame(data={"time": [], "lipid": [], "id": [], "resid_start": [], "resid_end": [], "leaflet": [],  "x": [], "y": [], "z": []})
            # convert column "time", "id", "resid_start", and "resid_end" to int64 dtype
            df_com = df_com.astype({"time": int, "id": int, "resid_start": int, "resid_end": int})
                
            filename = structpathname+"frame_dump_" + str(file_count) + ".gro"
            if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
                break
        
            tempora = round(initial_time + (file_count)*(vmdstep)*(dt)*(nstout))
            stdout.write("\n Taking care of snapshot %s for time %i ps \n" % (filename,tempora))
            # Find the box size of each gro file and append to boxsize.csv
            boxsize_read_gro(filename,continuationflag,tempora)
            
            # Find all bead locations for all lipids and store dataframes of positions in dictionary
            # Python Dictionary titled by LipidNames
            df_file = {}
            for n, name in enumerate(LipidNames):
                df_xyz = beads_read_gro(filename,name,beads_list[n])
                df_file[name] = df_xyz

            ############## CoM for each lipid in BackBoneBeads list CALCULATION ###################################
            stdout.write("Starting center of mass calculation for all specified beads in each lipid \n")

            for mol1,molecul1 in enumerate(LipidNames):
                df_xyz = df_file[molecul1]
                for lip in range(len(df_xyz)):
                    # vector sum of all beads designated for each lipid (row of dataframe)
                    vector = [0.0, 0.0, 0.0]
                    mass_sum = 0
                    for bbb in backbonebeads_list[mol1]:
                        vector[0] = vector[0] + float(df_beadmass[molecul1][bbb].iloc[0])*df_xyz.loc[lip, 'x_'+bbb]
                        vector[1] = vector[1] + float(df_beadmass[molecul1][bbb].iloc[0])*df_xyz.loc[lip, 'y_'+bbb]
                        vector[2] = vector[2] + float(df_beadmass[molecul1][bbb].iloc[0])*df_xyz.loc[lip, 'z_'+bbb]
                        mass_sum = mass_sum + float(df_beadmass[molecul1][bbb].iloc[0])
                    # center of mass of selected beads [assume beads are all weight: 72 amu...except CHOL]
                    # source for assumption: http://www.ks.uiuc.edu/Training/Tutorials/martini/rbcg-tutorial.pdf pg 11.
                    comvector = [x/mass_sum for x in vector]
                    df_com.loc[com_count, ['x','y','z']] = comvector
                    df_com.loc[com_count, ['time','lipid','id','resid_start', 'resid_end']] = [int(tempora), molecul1, int(df_xyz.loc[lip,'id']), int(df_xyz.loc[lip,'resid_start']), int(df_xyz.loc[lip,'resid_end'])]
                    com_count += 1
                # Save coordinate components of ATOC
                if molecul1 == 'ATOC':
                    datafilename = pathname+molecul1+'_beads_positions.dat'
                    saveordervalues(datafilename,df_xyz)
                # Save coordinate components of ATOC
                elif molecul1 == 'POPC':
                    datafilename = pathname+molecul1+'_beads_positions.dat'
                    saveordervalues(datafilename,df_xyz)
            ######################### Upper and Lower Leaflet Determination ########################################
            avg_z = df_com['z'].loc[df_com['time'] == tempora].mean()
            for lip in df_com['z'].loc[df_com['time'] == tempora].index:
                if df_com.loc[lip,'z'] >= avg_z:
                    df_com.loc[lip, 'leaflet'] = 'upper'
                elif df_com.loc[lip,'z'] < avg_z:
                    df_com.loc[lip,'leaflet'] = 'lower'
            ###################### INITIALIZE THETA, DX, DY, DZ PARAMETERS ########################################################
            for mol,molecul in enumerate(LipidNames):
                start_com = time.time()
                xydframe = df_com.loc[df_com['lipid'] == molecul] # need df_com info to build df_bond
                df_bond[molecul] = {}
                df_dx[molecul],df_dy[molecul],df_dz[molecul] = {},{},{}
                lipid_place = ""
                for i in xydframe['id']:
                    lipid_place+="%8i" % int(float(i))
                lipid_place+="\n"
                output_legend = "  Time[ps]" + lipid_place
                # Get list of individual bonds
                d = {'time': []}
                for pos in xydframe['id']:
                    d[str(int(float(pos)))] = []
                bonds = bond_list[mol].split()
                for bond in bonds:
                    df_bond[molecul][bond] = pd.DataFrame(data=d)   
                    df_dx[molecul][bond] = pd.DataFrame(data=d)
                    df_dy[molecul][bond] = pd.DataFrame(data=d)
                    df_dz[molecul][bond] = pd.DataFrame(data=d)
                #################################### SAVE COM TO TIMELOG #####################################################
                end_com = time.time()
                timelog('com',codefolderpath,molecul,mdnum,tempora, end_com, start_com)
                    
            ###################### (THETA CALCULATION #########################################################
            for mol,molecul in enumerate(LipidNames):
                start_theta = time.time()
                stdout.write("Starting theta calculation for all specified bonds in %s \n" % (molecul))
                # compute order parameter for each bond, for each snapshot
                df_xyz = df_file[molecul]
                bonds = bond_list[mol].split()
                for bond in df_bond[molecul]:
                    df_bond[molecul][bond].loc[file_count,'time'] = tempora
                    df_dx[molecul][bond].loc[file_count,'time'] = tempora
                    df_dy[molecul][bond].loc[file_count,'time'] = tempora
                    df_dz[molecul][bond].loc[file_count,'time'] = tempora
                for bo, bon in enumerate(bonds):
                    order_parameters = []
                    bond = bon.split('-')
                    first  = df_xyz.loc[:, ['x_'+bond[0], 'y_'+bond[0], 'z_'+bond[0]]]
                    second = df_xyz.loc[:, ['x_'+bond[1], 'y_'+bond[1], 'z_'+bond[1]]]
                    # compute order parameter for each lipid
                    for i in range(len(first)):
                        lipid_id = int(df_xyz.loc[i,'id'])
                        order_parameter = 0.0
                        # vector between the two previous beads (orientation doesn't matter)
                        vector = [0.0, 0.0, 0.0]
                        for j in range(3):
                            vector[j] = first.iloc[i,j] - second.iloc[i,j]
                        norm2 = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
                        # compute projection on the bilayer normal
                        projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
                        # order parameter
                        order_parameter = math.acos(projection/norm2)# convert strings to floats
                        # Append to dataframe in dictionary
                        df_bond[molecul][bon].loc[file_count,str(lipid_id)] = order_parameter
                        df_dx[molecul][bon].loc[file_count,str(lipid_id)] = vector[0]
                        df_dy[molecul][bon].loc[file_count,str(lipid_id)] = vector[1]
                        df_dz[molecul][bon].loc[file_count,str(lipid_id)] = vector[2]
                #################################### SAVE THETA TO TIMELOG #####################################################
                end_theta = time.time()
                timelog('theta',codefolderpath,molecul,mdnum,tempora, end_theta, start_theta)
            #################################### SAVE THETA, DX, DY, DZ #####################################################
            # SAVE TO FILES
            # write the center of mass results to csv file
            stdout.write("\n The momentaneous center of mass, theta, dx, dy, and dz of the selected lipids has been computed and saving. \n")
            saveordervalues(comstrFile, df_com)
            # write the theta, dx, dy, and dz to csv file
            for mol, molecul in enumerate(LipidNames):
                bonds = bond_list[mol].split()
                for bond in bonds:
                    datafilename = pathname+molecul+'_theta_'+bond+'.dat'
                    saveordervalues(datafilename,df_bond[molecul][bond])
                    datafilename = pathname+molecul+'_dx_'+bond+'.dat'
                    saveordervalues(datafilename,df_dx[molecul][bond])
                    datafilename = pathname+molecul+'_dy_'+bond+'.dat'
                    saveordervalues(datafilename,df_dy[molecul][bond])
                    datafilename = pathname+molecul+'_dz_'+bond+'.dat'
                    saveordervalues(datafilename,df_dz[molecul][bond])
                    
            # Remove .gro structure file (requires too much memory)
            remove(filename)
            
            # Next gro file
            file_count += 1
            

        # write the center of mass results to csv file
        stdout.write("\n All the center of mass of the selected lipids has been computed. \n\n")
       

        end_total = time.time()
        stdout.write("Runtime of CoM and theta calculation is %0.4f seconds \n" % (end_total - start_total))
        #################################### SAVE COM+THETA TO TIMELOG #####################################################
        timelog('comtotal',codefolderpath,'',mdnum,'', end_total, start_total)
        
        ############################### IMPORT THE FULL df_com  FILE to return df_com ############################# 
        start_com = time.time()
        import numpy as np
        file1 = open(comstrFile,'r')
        Lines = file1.readlines()
        # separate spaces between string data points
        strLines = [Lines[i].split(',') for i in range(0,len(Lines))]
        # remove '\n' from lines
        for i, sline in enumerate(strLines):
            for j, svalue in enumerate(sline):
                strLines[i][j] = svalue.strip('\n')
        # convert into arrays
        strArray = np.array(strLines)
        # create dataframe from array with columns defined by title line
        df_com = pd.DataFrame(strArray[1:], columns=strArray[0])
        # convert strings to floats
        df_com = df_com.astype({'id': 'float', 'resid_start': 'float', 'resid_end': 'float', 'time': 'float','x':'float','y':'float','z':'float'})
        # convert floats to integers
        df_com = df_com.astype({'id': 'int', 'resid_start': 'int', 'resid_end': 'int', 'time': 'int'})
        end_com = time.time()
        stdout.write("CoM import runtime is %0.4f seconds \n\n" % (end_com - start_com))
        
    else:
        ############################### import the already calculated df_com  file to return df_com ############################# 
        start_com = time.time()
        import numpy as np
        stdout.write("The center of mass file is: " + comstrFile+'\n')
        stdout.write("Continuing with center of mass import... \n")
        #df_com =pd.read_csv(comstrFile,index_col=False) 
        file1 = open(comstrFile,'r')
        Lines = file1.readlines()

        # separate spaces between string data points
        strLines = [Lines[i].split(',') for i in range(0,len(Lines))]
        # remove '\n' from lines
        for i, sline in enumerate(strLines):
            for j, svalue in enumerate(sline):
                strLines[i][j] = svalue.strip('\n')
        # convert into arrays
        strArray = np.array(strLines)
        # create dataframe from array with columns defined by title line
        df_com = pd.DataFrame(strArray[1:], columns=strArray[0])
        # convert strings to floats
        df_com = df_com.astype({'id': 'float', 'resid_start': 'float', 'resid_end': 'float', 'time': 'float','x':'float','y':'float','z':'float'})
        # convert floats to integers
        df_com = df_com.astype({'id': 'int', 'resid_start': 'int', 'resid_end': 'int', 'time': 'int'})
        end_com = time.time()
        stdout.write("CoM import runtime is %0.4f seconds \n\n" % (end_com - start_com))
    return df_com
