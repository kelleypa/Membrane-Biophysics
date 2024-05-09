#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2022.02.23 V7 generalize for all lipids & combine all analysis - Patrick Kelley


##########################################################################################################
################# PBC FLOOD FILL FOR DOMAIN DETERMINATION  ###############################################
##########################################################################################################

# NEAREST 1 NEIGHBOR (DIAGONAL) FLOOD FILL ALGORITHM
def fillrin(matrix, domaintype, x, y, advnum, area):
    import numpy as np
    def fillrin(x, y, area):
        if x < 0 or x >= matWidth or y < 0 or y >= matHeight:
            if x < 0:
                area = fillrin(matHeight-1,y,area)
            if y < 0:
                area = fillrin(x,matWidth-1,area)
            if x >= matWidth:
                area = fillrin(0,y,area)
            if y >= matHeight:
                area = fillrin(x,0,area)
        elif 0 <= x < matWidth and 0 <= y < matHeight:
            if matrix[y][x] == domaintype:
                matrix[y][x] = advnum
                area += 1

                area = fillrin(x, y - 1, area) #south
                area = fillrin(x, y + 1, area) #north
                area = fillrin(x + 1, y, area) #west
                area = fillrin(x - 1, y, area) #east
                area = fillrin(x + 1, y - 1, area) #southwest
                area = fillrin(x - 1, y + 1, area) #northeast
                area = fillrin(x + 1, y + 1, area) #northwest
                area = fillrin(x - 1, y - 1, area) #southeast
              
        return area
    matWidth = np.size(matrix,1)
    matHeight = np.size(matrix,0)
    area = fillrin(x, y, area)
    return matrix,area

# NEAREST 2 NEIGHBOR (DIAGONAL) FLOOD FILL ALGORITHM
def fill2rin(matrix, domaintype, x, y, advnum, area):
    import numpy as np
    def fill2rin(x, y, area):
        if x < 0 or x >= matWidth or y < 0 or y >= matHeight:
            if x < 0:
                area = fill2rin(matHeight-1,y,area)
            if y < 0:
                area = fill2rin(x,matWidth-1,area)
            if x >= matWidth:
                area = fill2rin(0,y,area)
            if y >= matHeight:
                area = fill2rin(x,0,area)
        elif 0 <= x < matWidth and 0 <= y < matHeight:
            if matrix[y][x] == domaintype:
                matrix[y][x] = advnum
                area += 1

                area = fill2rin(x, y - 1, area) #south
                area = fill2rin(x, y + 1, area) #north
                area = fill2rin(x + 1, y, area) #west
                area = fill2rin(x - 1, y, area) #east
                area = fill2rin(x + 1, y - 1, area) #southwest
                area = fill2rin(x - 1, y + 1, area) #northeast
                area = fill2rin(x + 1, y + 1, area) #northwest
                area = fill2rin(x - 1, y - 1, area) #southeast
                # Nearest second neighbor
                area = fill2rin(x, y - 2, area)
                area = fill2rin(x, y + 2, area)
                area = fill2rin(x - 1, y - 2, area)
                area = fill2rin(x - 1, y + 2, area)
                area = fill2rin(x + 1, y - 2, area)
                area = fill2rin(x + 1, y + 2, area)
                area = fill2rin(x - 2, y + 2, area)
                area = fill2rin(x - 2, y + 1, area)
                area = fill2rin(x - 2, y - 1, area)
                area = fill2rin(x - 2, y, area)
                area = fill2rin(x - 2, y - 2, area)
                area = fill2rin(x + 2, y - 2, area)
                area = fill2rin(x + 2, y - 1, area)
                area = fill2rin(x + 2, y, area)
                area = fill2rin(x + 2, y + 1, area)
                area = fill2rin(x + 2, y + 2, area)
        return area
    matWidth = np.size(matrix,1)
    matHeight = np.size(matrix,0)
    area = fill2rin(x, y, area)
    return matrix,area

def threebythreebox(points, x_box, y_box):
    import pandas as pd
    point00 = points.copy()
    point10 = points.copy()
    point20 = points.copy()
    point01 = points.copy()
    point02 = points.copy()
    point11 = points.copy()
    point12 = points.copy()
    point21 = points.copy()
    point22 = points.copy()

    point00['x'] = point00['x'] - x_box
    point00['y'] = point00['y'] - y_box
    point10['x'] = point10['x'] - x_box
    point20['x'] = point20['x'] - x_box
    point20['y'] = point20['y'] + y_box
    point01['y'] = point01['y'] - y_box
    point02['x'] = point02['x'] + x_box
    point02['y'] = point02['y'] - y_box
    point12['x'] = point12['x'] + x_box
    point21['y'] = point21['y'] + y_box
    point22['x'] = point22['x'] + x_box
    point22['y'] = point22['y'] + y_box
    
    threebythreebox = pd.concat([point00,point10,point20,point01,point02,point11,point12,point21,point22])
    
    return threebythreebox
def find_nearest(a, a0):
    import numpy as np
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(np.array(a) - a0).argmin()
    return idx, a.flat[idx]
def find_csv_filenames( path_to_dir, suffix=".csv" ):
    from os import listdir
    filenames = listdir(path_to_dir)
    return [ filename for filename in filenames if filename.endswith( suffix ) ]
def beads_read_gro(file, molecule, atoms):
    import pandas as pd
    df_name = []
    df_name.append('id')
    df_name.append('resid_start')
    df_name.append('resid_end')
    for atomic in atoms:
        df_name.append('x_'+atomic)
        df_name.append('y_'+atomic)
        df_name.append('z_'+atomic)
    d = {pos: [] for pos in df_name}
    df_xyz = pd.DataFrame(data = d)

    number_of_particles = 0
    lipid_bead_count = 0
    lipid_row = 0
    line_counter = 0
    for line in open(file):
        if line_counter == 1:
            number_of_particles = int(line)
        elif line_counter > 1 and line_counter < number_of_particles + 2:
            atomic = line[10:15].strip()
            molecular = line[5:10].strip()
            if atomic in atoms and molecular in molecule: 
                # id and starting residue number
                if lipid_bead_count == 0: 
                    df_xyz.loc[lipid_row,'id'] = float(line[0:5].strip()) # id number
                    df_xyz.loc[lipid_row,'resid_start'] = float(line[15:20].strip()) # start residue number
                # bead positions
                df_xyz.loc[lipid_row,'x_'+atomic] = float(line[20:28])
                df_xyz.loc[lipid_row,'y_'+atomic] = float(line[28:36])
                df_xyz.loc[lipid_row,'z_'+atomic] = float(line[36:44])
                # update number of bead in Molecule
                lipid_bead_count += 1
                # last bead and last residue number
                if lipid_bead_count == len(atoms):
                    df_xyz.loc[lipid_row,'resid_end'] = float(line[15:20].strip()) # end residue number
                    lipid_row += 1
                    lipid_bead_count = 0
        line_counter += 1
    return df_xyz
def boxsize_read_gro(filename,continuationflag,tempora):
    import pandas as pd
    import os
    from sys import argv, stdout
    def isfloat(num):
        try:
            float(num)
            return True
        except ValueError:
            return False
    boxpath = "boxsize.csv"
    # open the file
    f = open(filename, "r")
    # read all the lines
    lines = f.readlines() 
    # close the file
    f.close()

    # TIME: get the time in the first row of the gro file
    if continuationflag == True:
        tst = [int(float(i)) for i in lines[0].split() if isfloat(i)]
    else:
        tst = [tempora]

    # BOXSIZE: split last line, ie the box vectors, in gro file by empty spaces
    lst = lines[-1].strip().split(' ')
    # remove empty elements in list
    lst = list(filter(None, lst))
    # convert str to float
    lst = [float(i) for i in lst] 

    # create dataframe of gro files boxsize
    df_box = pd.DataFrame(data={'t[ps]':tst[0],'x':lst[0],'y':lst[1],'z':[lst[2]]})
    # if boxsize.csv does exist, append new time boxsize elements to last known entry
    if os.path.exists(boxpath) == True:
        # read boxsize data already in boxpath
        df_boximport = pd.read_csv(boxpath)
        # if new continuous time is already computed
        if tst[0] == df_boximport['t[ps]'].iloc[-1]:
            stdout.write('time already exists, continuing \n')
        elif tst[0] > df_boximport['t[ps]'].iloc[-1]:
            with open(boxpath, mode='a', newline='') as boxcsvfile:
                df_box.to_csv(boxcsvfile, header=(boxcsvfile.tell()==0), index=False)
        else:
            # merge imported boxsize data with gro file boxsize
            df_merged = pd.concat([df_boximport,df_box],axis=0)
            # sort and remove duplicates by time
            df_boxsize = df_merged.sort_values(by=['t[ps]']).drop_duplicates('t[ps]', keep='first')
            # save new to csv
            with open(boxpath, mode='w', newline='') as boxcsvfile:
                df_boxsize.to_csv(boxcsvfile, header=(boxcsvfile.tell()==0), index=False)
    # if boxsize.csv does exist but time is not continuous, append elements without checking last known entry
    # if boxsize.csv doesn't exist, create it and add first element
    else:
        with open(boxpath, mode='a', newline='') as boxcsvfile:
            df_box.to_csv(boxcsvfile, header=(boxcsvfile.tell()==0), index=False)

def savefigurestimesincrements(tempura,figsaveinterval):
    import numpy as np
    def find_nearest_time(a, a0):
        import numpy as np
        "Element in nd array `a` closest to the scalar value `a0`"
        idx = np.abs(np.asarray(a) - a0).argmin()
        return a[idx]
    figsavevery = np.arange(tempura[0],tempura[-1],figsaveinterval)
    saveveryfig = []
    for figgy in figsavevery:
        saveveryfig.append(find_nearest_time(tempura, figgy))
    return saveveryfig
##########################################################################################################
#################### Variable Threshold Values ###########################################################
##########################################################################################################
def variablethresholdvalues(maxwindowlipid, lipidchar,lipidnumber,randperc,raftlike,nonraftlike,codefolderpath):
    import math
    import numpy as np
    from itertools import combinations
    from itertools import combinations_with_replacement
    import os
    import csv
    import pandas as pd
    from sys import stdout
    ############## BINOMIAL COEFFICIENTS ##################################
    def binom(N,m):
        binomcoeff = math.factorial(N)/(math.factorial(m)*math.factorial(N-m))
        return binomcoeff
    #######################################################################
    totlip = np.sum(lipidnumber)
    raftdenscnt,nonraftdenscnt = [], []
    raftcholdpsmpercthres, nonraftphospercthres = [], []
    
    # Run from last variable thresholdvalue
    vtvlogpath = codefolderpath+'variablethresholdvalues.csv'
    isExist = os.path.exists(vtvlogpath)
    if isExist:
        # append data frame to CSV file
        df = pd.read_csv(vtvlogpath, sep=',', index_col=None)
        last_line = df.rlist.values[-1]
        if last_line == maxwindowlipid:
            stdout.write("\n\nAll " + str(maxwindowlipid)+ " Variable Threshold Values had been computed. Continuing...\n\n")
            raftcholdpsmpercthres = df.raftpercthres
            nonraftphospercthres = df.nonraftpercthres
            return raftcholdpsmpercthres, nonraftphospercthres
        else:
            stdout.write("\nVariable Threshold Value up to " + str(last_line)+" out of " + str(maxwindowlipid) + ". Computing further.\n\n")
            rlist = range(last_line+1,maxwindowlipid+1,1)
    elif not isExist: 
        rlist = range(1,maxwindowlipid+1,1)
    
    for r in rlist:
        ###### Create Combination with Repetition List ######
        combo = list(combinations_with_replacement(lipidchar, r))
        ### Find Weighted Probability of each Microstate ####
        probweight,CDcntlist,PLcntlist = [],[],[]

        for comb in combo:
            countlist = [0]*len(lipidnumber)
            indexlist = []
            CDcnt = 0
            PLcnt = 0
            for i in comb:
                index = lipidchar.index(i)
                indexlist = indexlist + [index]
                countlist[index]+=1

                if i in raftlike:
                    CDcnt+=1
                elif i in nonraftlike:
                    PLcnt+=1
            CDcntlist.append(CDcnt/r)
            PLcntlist.append(PLcnt/r)
            probw = 1
            for c, cnt in enumerate(countlist):
                probw = probw * binom(lipidnumber[c],cnt)
            probweight = probweight + [probw/binom(totlip,r)]

        # using list comprehension to get numbers > k
        phospercthres = False
        chsmpercthres = False
        startpercthres = True
        for percthres in np.arange(0,1,0.001)[::-1]:
            if chsmpercthres == True and phospercthres == True:
                break
            elif startpercthres == True:
                Pcountprev  = [probweight[i] for i,PL in enumerate(PLcntlist) if float(PL) >= float(percthres)]
                CDcountprev = [probweight[i] for i,CD in enumerate(CDcntlist) if float(CD) >= float(percthres)]
                percthresprev = percthres
                startpercthres = False
            else:
                Pcount  = [probweight[i] for i,PL in enumerate(PLcntlist) if float(PL) >= float(percthres)]
                CDcount = [probweight[i] for i,CD in enumerate(CDcntlist) if float(CD) >= float(percthres)]
                # Nonraft (Phospholipid)
                if np.sum(Pcount) > (randperc) and phospercthres == False: # once less than percentage chance
                    nonraftdenscnt.append(round(np.sum(Pcountprev),4))
                    nonraftphospercthres.append(percthresprev)
                    phospercthres = True
                # Raft (Cholesterol/ DPSM)
                if (np.sum(CDcount)) > (randperc) and chsmpercthres == False: # once less than percentage chance
                    raftdenscnt.append(round(np.sum(CDcountprev),4))
                    raftcholdpsmpercthres.append(percthresprev)
                    chsmpercthres = True
                Pcountprev  = Pcount
                CDcountprev = CDcount
                percthresprev = percthres

    # dictionary of lists 
    dict = {'rlist': rlist, 'raftpercthres': raftcholdpsmpercthres, 'nonraftpercthres': nonraftphospercthres} 
    df = pd.DataFrame(dict)
    # append data frame to CSV file
    hdr = False  if os.path.isfile(vtvlogpath) else True
    df.to_csv(vtvlogpath, mode='a', index=False, header=hdr)
    
    return raftcholdpsmpercthres, nonraftphospercthres

#########################################################################################################
############################################ PLOTS ######################################################
#########################################################################################################
def leafletdomainplot(leaf,picpath,phospholipid,LipidNames,temp,randperc,figsave,leafxi,leafyi,leafboxintplot,leafpoints,leafspeed,pointlnwth,pointsize,x_box,y_box,xwindow,ywindow,xstep,ystep):

    #################################### PLOT SETTINGS ######################################################
    #MATPLOTLIB 
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.ticker as mtick
    import matplotlib.cm as cm
    import matplotlib.patches as mpatches
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib.patches import Patch
    from matplotlib.legend_handler import HandlerPatch
    # SOME COLORS
    purple = "#7E1E9C"
    violet = "#9A0EEA"
    darkblue = "#030764"

    class AnyObject(object):
        pass
    class AnyObjectHandler(object):
        def create_artists(self, legend, orig_handle, fontsize, handlebox):
            x0, y0 = handlebox.xdescent, handlebox.ydescent
            width, height = handlebox.width, handlebox.height
            p = mpatches.Rectangle([x0, y0], width, height,
                                    transform=handlebox.get_transform())
            handlebox.add_artist(p)
            return p
    plt.style.use('default')
    #########################################################################################################
    pltfilename = picpath+phospholipid+'_domain_'+leaf+'_'+str(int(temp))+'ps.png'
    import os
    isExist = os.path.exists(pltfilename)
    if not isExist and figsave == True:
        # save and clear
        fig, ax1 = plt.subplots( figsize=(14,12))
        raftypoints = False
        ##########################################################################################################
        # LEAFLET CONTOUR PLOT
        cmap = colors.ListedColormap(['cornflowerblue', 'white', 'tomato'])
        CS1 = ax1.contourf(leafxi, leafyi, leafboxintplot, 2, cmap=cmap ,
                        vmax=2, vmin=0)

        # LEAFLET SCATTER PLOT
        ax1.scatter(leafpoints[:, 0], leafpoints[:, 1],c = leafspeed, edgecolors='black', linewidth=pointlnwth, s=pointsize)
        ax1.set_aspect("equal")
        ax1.set_title(('Voting: '+'t='+str(int(temp))+'ps, '+r'Variable Density Threshold (P $\geq$ %0.1f%%)'+'\n'+r'Window = %0.1f$\AA$ $\times$ %0.1f$\AA$, $\Delta x_{slide}$=%0.2f$\AA$, $\Delta y_{slide}$=%0.2f $\AA$'+'\n'+r'Leaflet: %s, Box = %0.1f$\AA$$\times$%0.1f$\AA$')%(randperc*100,xwindow*10,ywindow*10,xstep*10,ystep*10,leaf,x_box*10,y_box*10),fontsize=25)

        # COLORBAR
        cbar1 = fig.colorbar(CS1,ax=[ax1], pad = 0.01)
        cbar1.set_ticks(range(0,3))
        cbar1.set_ticklabels(['PL-rich', 'Mixed', 'CHOL/SM-rich'])
        # X & Y LABEL AND FONTSIZE
        ax1.set_xlabel('x [nm]', fontsize = 20)
        ax1.set_ylabel('y [nm]', fontsize = 20, rotation=90, labelpad=10)
        ax1.FontSize = 25; 

        # LEGEND
        if raftypoints == False:
            patchy = []
            for lip in LipidNames:
                if lip == 'DPSM':
                    p = Patch(facecolor='r',
                            edgecolor='black', linewidth=3)
                elif lip == 'CHOL':
                    p = Patch(facecolor='w',
                            edgecolor='black', linewidth=3)
                elif lip == 'POPC':
                    p = Patch(facecolor='b',
                            edgecolor='black', linewidth=3)
#                     elif lip == 'PUPC':
#                         p = Patch(facecolor=darkblue,
#                                 edgecolor='black', linewidth=3)
                elif lip == 'PUPC':
                    p = Patch(facecolor='blue',
                            edgecolor='black', linewidth=3)
              
# Colors: Orange - 'orange', Amber - '#FFBF00', Yellow - 'yellow', Pink - 'pink', Green - 'green', Black - 'black'
                elif lip == 'ATOC':
                    p = Patch(facecolor='#FFBF00', 
                            edgecolor='black', linewidth=3)
                patchy = patchy+[p]
            fig.legend(patchy, LipidNames, handler_map={AnyObject: AnyObjectHandler()},prop={'size': 20},bbox_to_anchor=(0.25,-0.01,0.5,0.5), loc="lower left", mode="expand", borderaxespad=0, ncol=len(LipidNames))

        plt.savefig(pltfilename)
        plt.close()
def leafletscatterplot(leaf,picpath,phospholipid,LipidNames,temp,randperc,figsave,leafpoints,leafspeed,pointlnwth,pointsize,x_box,y_box,xwindow,ywindow,xstep,ystep):
    
    #################################### PLOT SETTINGS ######################################################
    #MATPLOTLIB 
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.ticker as mtick
    import matplotlib.cm as cm
    import matplotlib.patches as mpatches
    from matplotlib_scalebar.scalebar import ScaleBar
    from matplotlib.patches import Patch
    from matplotlib.legend_handler import HandlerPatch
    # SOME COLORS
    purple = "#7E1E9C"
    violet = "#9A0EEA"
    darkblue = "#030764"

    class AnyObject(object):
        pass
    class AnyObjectHandler(object):
        def create_artists(self, legend, orig_handle, fontsize, handlebox):
            x0, y0 = handlebox.xdescent, handlebox.ydescent
            width, height = handlebox.width, handlebox.height
            p = mpatches.Rectangle([x0, y0], width, height,
                                    transform=handlebox.get_transform())
            handlebox.add_artist(p)
            return p
    plt.style.use('default')
    #########################################################################################################
    pltfilename = picpath+phospholipid+'_scatter_'+leaf+'_'+str(int(temp))+'ps.png'
    import os
    isExist = os.path.exists(pltfilename)
    if not isExist and figsave == True:
        # save and clear
        fig, ax1 = plt.subplots( figsize=(14,12))
        raftypoints = False
        ##########################################################################################################

        # LEAFLET SCATTER PLOT  
        ax1.scatter(leafpoints[:, 0], leafpoints[:, 1],c = leafspeed,edgecolors='black', linewidth=pointlnwth, s=pointsize)
        ax1.set_facecolor("lightgrey")

        ax1.set_aspect("equal")
        ax1.set_title(('Voting: '+'t='+str(int(temp))+'ps, '+r'Variable Density Threshold (P $\geq$ %0.1f%%)'+'\n'+r'Window = %0.1f$\AA$ $\times$ %0.1f$\AA$, $\Delta x_{slide}$=%0.2f$\AA$, $\Delta y_{slide}$=%0.2f $\AA$'+'\n'+r'Leaflet: %s, Box = %0.1f$\AA$$\times$%0.1f$\AA$')%(randperc*100,xwindow*10,ywindow*10,xstep*10,ystep*10,leaf,x_box*10,y_box*10),fontsize=25)

        # X & Y LABEL AND FONTSIZE
        ax1.set_xlabel('x [nm]', fontsize = 20)
        ax1.set_ylabel('y [nm]', fontsize = 20, rotation=90, labelpad=10)
        ax1.FontSize = 25; 

        # LEGEND
        if raftypoints == False:
            patchy = []
            for lip in LipidNames:
                if lip == 'DPSM':
                    p = Patch(facecolor='r',
                            edgecolor='black', linewidth=3)
                elif lip == 'CHOL':
                    p = Patch(facecolor='w',
                            edgecolor='black', linewidth=3)
                elif lip == 'POPC':
                    p = Patch(facecolor='b',
                            edgecolor='black', linewidth=3)
#                     elif lip == 'PUPC':
#                         p = Patch(facecolor=blue,
#                                 edgecolor='black', linewidth=3)
                elif lip == 'PUPC':
                    p = Patch(facecolor='darkblue',
                            edgecolor='black', linewidth=3)
                elif lip == 'ATOC':
                    p = Patch(facecolor='#FFBF00',
                            edgecolor='black', linewidth=3)
                patchy = patchy+[p]
            fig.legend(patchy, LipidNames, handler_map={AnyObject: AnyObjectHandler()},prop={'size': 20},bbox_to_anchor=(0.25,-0.01,0.5,0.5), loc="lower left", mode="expand", borderaxespad=0, ncol=len(LipidNames))

        plt.savefig(pltfilename)
        plt.close()

##########################################################################################################
# TIME LOG
##########################################################################################################
def timelog(timelogcomp,codefolderpath,phospholipid,mdnum,temp, end_time, start_time):
    timelogpath = codefolderpath+phospholipid+'_timelog.txt'
    
    if timelogcomp == 'com': #center of mass
        with open(timelogpath, "a") as f:
            output_legend = 80
            f.write(("-"*(output_legend)) + "\n")
            f.write("Runtime of %s CoM Calculation for MD%i at %i ps: %0.4f seconds \n" % (phospholipid,mdnum,temp, end_time - start_time))
            f.close()
    elif timelogcomp == 'cos2': # domain determination
        with open(timelogpath, "a") as f:
            output_legend = 80
            f.write(("-"*(output_legend)) + "\n")
            f.write("Runtime of %s cos2 for MD%i at %i ps: %0.4f seconds \n" % (phospholipid,mdnum,temp, end_time - start_time))
            f.close()
    elif timelogcomp == 'comtotal': # domain determination
        with open(timelogpath, "a") as f:
            output_legend = 90
            f.write(("-"*(output_legend)) + "\n")
            f.write(("-"*(output_legend)) + "\n")
            f.write("Total Runtime of CoM and cos2 Calculation for MD%i: %0.4f seconds \n" % (mdnum, end_time - start_time))
            f.write(("-"*(output_legend)) + "\n")
            f.close()
    elif timelogcomp == 'dd': # domain determination
        with open(timelogpath, "a") as f:
            output_legend = 80
            f.write(("-"*(output_legend)) + "\n")
            f.write("Runtime of %s Domain Determination for MD%i at %i ps: %0.4f seconds \n" % (phospholipid,mdnum,temp, end_time - start_time))
            f.close()
##########################################################################################################
##########################################################################################################
# CENTER OF MASS / ORDER PARAMETER FILE SAVING
##########################################################################################################
##########################################################################################################
def saveordervalues(datapath,df):
    import os
    import csv
    isExist = os.path.exists(datapath)
    if isExist:
        # append data frame to CSV file
        df.to_csv(datapath, sep=',', mode='a', index=False, header=False)
    elif not isExist:
        # create the new data file
        with open(datapath, 'w') as csvfile:
            # using csv.writer method from CSV package
            write = csv.writer(csvfile)
            write.writerow(df.columns.values)
            # Close the file object
            csvfile.close()
        df.to_csv(datapath, sep=',', mode='a', index=False, header=False)

##########################################################################################################
##########################################################################################################
# DOMAIN DETERMINATION MAIN ALGORITHM
##########################################################################################################
##########################################################################################################
def leafletdomaindetermine(df_com,leaf,temp, xstep,ystep,xwindow,ywindow,xsteplist,ysteplist, x_box,y_box, raftcholdpsmpercthres, nonraftphospercthres):
    import pandas as pd
    import numpy as np
    import math, copy
    from functions import fillrin, fill2rin, threebythreebox, find_nearest
    # plots
    import itertools
    from scipy.interpolate import griddata
    import sys
    sys.setrecursionlimit(6000)
    from sys import argv, stdout
    
    df_temp = df_com.loc[(df_com['time'] == temp) & (df_com['leaflet'] == leaf)].copy()
    
    gridbox_y = np.arange(ystep/2,y_box+ystep/2,ystep)
    gridbox_x = np.arange(xstep/2,x_box+xstep/2,xstep)

    # for Nearest Neighbor Interpolation
    denX,denY,denZ = [],[],[]
    domain = np.empty([len(np.arange(ystep/2,y_box+ystep/2,ystep)),len(np.arange(xstep/2,x_box+xstep/2,xstep))])
    domain[:,:] = np.nan
    # for Voting 
    raftvote = copy.deepcopy(domain)
    raftvote[:,:] = 0
    nonraftvote = copy.deepcopy(raftvote)
    mixedvote = copy.deepcopy(raftvote)

    df_temp3x3 = threebythreebox(copy.deepcopy(df_temp),x_box, y_box)
    
    for iy, dy in enumerate(ysteplist):
        dy1 = round(dy-ywindow+ystep,1)
        dy2 = round(dy+ystep,1)
        df_y = df_temp3x3.loc[(df_temp3x3['y'] >= dy1) & (df_temp3x3['y'] <= dy2)].copy()
        df_y = df_y.sort_values(by=['x'])
        df_choldpsm = df_y.loc[(df_y['lipid'] == 'CHOL') | (df_y['lipid'] == 'DPSM')].copy()
        df_pl = df_y.loc[(df_y['lipid'] == 'POPC') | (df_y['lipid'] == 'PUPC')].copy()
        for ix, dx in enumerate(xsteplist):
            dx1 = round(dx-xwindow+xstep,1)
            dx2 = round(dx+xstep,1)
            df_choldpsmcount = df_choldpsm.loc[(df_choldpsm['x'] >= dx1) & (df_choldpsm['x'] <= dx2)]
            df_plcount = df_pl.loc[(df_pl['x'] >= dx1) & (df_pl['x'] <= dx2)]
            df_x = df_y.loc[(df_y['x'] >= dx1) & (df_y['x'] <= dx2)]
            choldpsmcount = df_choldpsmcount['lipid'].count()
            plcount = df_plcount['lipid'].count()
            windowcount = df_x['lipid'].count()
            # for Nearest Neighbor Interpolation
            denX.append(dx)
            denY.append(dy)

            raftrashio = choldpsmcount/windowcount
            nonraftrashio = plcount/windowcount
            # for Nearest Neighbor Interpolation
            denZ.append(raftrashio)

            indy = np.where((gridbox_y >= dy1) & (gridbox_y <= dy2))[0]
            indx = np.where((gridbox_x >= dx1) & (gridbox_x <= dx2))[0]
            # TESTING
            # stdout.write(str(len(indy)) + '\t' + str(len(indx)) + '\t' + str(windowcount) +'\n')
            # stdout.write(str(indy) + '\t' + str(indx) +'\n')
            if windowcount == 0:
                continue
            elif raftrashio > raftcholdpsmpercthres[windowcount-1]: #raft with variable threshold ratio
                for vy in range(indy[0],indy[-1]+1,1):
                    for vx in range(indx[0],indx[-1]+1,1):
                        raftvote[vy,vx] += 1
            elif nonraftrashio > nonraftphospercthres[windowcount-1]:
                for vy in range(indy[0],indy[-1]+1,1):
                    for vx in range(indx[0],indx[-1]+1,1):
                        nonraftvote[vy,vx] += 1
            else:
                for vy in range(indy[0],indy[-1]+1,1):
                    for vx in range(indx[0],indx[-1]+1,1):
                        mixedvote[vy,vx] += 1

    ##########################################################################################################
    # LEAFLET: VOTING RAFT, NONRAFT, INTERMEDIATE  
    ##########################################################################################################
    # https://towardsdatascience.com/8-ways-to-filter-pandas-dataframes-d34ba585c1b8?gi=ca9fe6dcdf32
    boxint = copy.deepcopy(domain)
    for iy in range(0,np.size(domain,0)):
        for ix in range(0,np.size(domain,1)):
            totalvotecnt = raftvote[iy,ix]+nonraftvote[iy,ix]+mixedvote[iy,ix]
            # TESTING
            #stdout.write(str(raftvote[iy,ix]) + '\t'+ str(nonraftvote[iy,ix]) + '\t'+ str(mixedvote[iy,ix]) + '\t' + str(int(totalvotecnt)) + '\n')
            raftvotecnt = raftvote[iy,ix]/totalvotecnt
            nonraftvotecnt = nonraftvote[iy,ix]/totalvotecnt
            mixedvotecnt = mixedvote[iy,ix]/totalvotecnt
            if raftvotecnt > nonraftvotecnt and raftvotecnt > mixedvotecnt: #raft
                boxint[iy,ix] = 2
            elif nonraftvotecnt > raftvotecnt and nonraftvotecnt > mixedvotecnt: #nonraft
                boxint[iy,ix] = 0
            else: #mixed
                boxint[iy,ix] = 1
    ####################################################################################
    #################### FIND LIPIDS IN EACH GRID ######################################
    ####################################################################################
    df_temp = df_com.loc[(df_com['time'] == temp) & (df_com['leaflet'] == leaf)].copy()
    ids = []
    ygridsteplist = np.arange(0,y_box+ystep,ystep)
    xgridsteplist = np.arange(0,x_box+xstep,xstep)
    for iy, dy in enumerate(ygridsteplist[0:-1]):
        dy1 = round(dy,3)
        dy2 = round(ygridsteplist[iy+1],3)
        # lipids outside of bounds of box
        if dy1 == round(0,3):
            df_y = df_temp.loc[(df_temp['y'] <= dy2)].copy()
        elif dy2 == round(ygridsteplist[-1],3):
            df_y = df_temp.loc[(df_temp['y'] >= dy1)].copy()
        else:
            df_y = df_temp.loc[(df_temp['y'] >= dy1) & (df_temp['y'] <= dy2)].copy()
        df_y = df_y.sort_values(by=['x'])
        for ix, dx in enumerate(xgridsteplist[0:-1]):
            dx1 = round(dx,3)
            dx2 = round(xgridsteplist[ix+1],3)
            # lipids outside of bounds of box
            if dx1 == round(0,3):
                df_x = df_y.loc[(df_y['x'] <= dx2)]
            elif dx2 == round(xgridsteplist[-1],3):
                df_x = df_y.loc[(df_y['x'] >= dx1)]
            else:
                df_x = df_y.loc[(df_y['x'] >= dx1) & (df_y['x'] <= dx2)]

            ids.append(df_x['id'].values.tolist())
    ####################################################################################
    ############## IDENTIFY LIPIDS IN EACH DOMAIN ######################################
    ####################################################################################        
    raftids,intermedids,nonraftids = [],[],[]
    advnum = 0
    for iy, dy in enumerate(boxint[:,0]):
        for ix, dx in enumerate(boxint[0,:]):
            if boxint[iy][ix] == 0:
                nonraftids.append(ids[advnum])
            elif boxint[iy][ix] == 1:
                intermedids.append(ids[advnum])
            elif boxint[iy][ix] == 2:
                raftids.append(ids[advnum])
            advnum += 1
    # Flatten raft list of lists
    flattened_raftids = list(itertools.chain(*raftids))
    raftidslist = copy.deepcopy(raftids)
    raftids = [*set(flattened_raftids)]
    # Flatten nonraft list of lists
    flattened_nonraftids = list(itertools.chain(*nonraftids))
    nonraftidslist = copy.deepcopy(nonraftids)
    nonraftids = [*set(flattened_nonraftids)]
    # Flatten intermediate list of lists
    flattened_intermedids = list(itertools.chain(*intermedids))
    intermedidslist = copy.deepcopy(intermedids)
    intermedids = [*set(flattened_intermedids)]

    raftpoints = df_temp[df_temp.id.isin(raftids)]
    nonraftpoints = df_temp[df_temp.id.isin(nonraftids)]
    interpoints = df_temp[df_temp.id.isin(intermedids)]

    xydframe = df_com.loc[df_com['time'] == temp]
    xydframe = xydframe.loc[xydframe['leaflet'] == leaf] # coordinates of all lipids in select frame
    xydframe.index = range(len(xydframe)) # set the indices to 0 of upper or lower leaflet 
    xydframe = xydframe.loc[:,['lipid','id','x','y','z']] # x, y, z coordinate for thickness map
    points = xydframe.loc[:,['x','y']].to_numpy()
    ##########################################################################################################
    # Convert from pandas dataframes to numpy arrays
    denX, denY, denZ = np.array(denX), np.array(denY), np.array(denZ)
    ################################# BLOCKY INTERPOLATION ###################################################
    # create x-y points to be used in heatmap
    # if static grid (sliding full length of window); else sliding fraction of window
    if xstep == xwindow and ystep==ywindow:
        xi = np.linspace(denX.min(), np.arange(0,x_box+xstep,xstep)[-1], 10*np.size(domain,1))
        yi = np.linspace(denY.min(), np.arange(0,y_box+ystep,ystep)[-1], 10*np.size(domain,0))
    else:    
        xi = np.linspace(denX.min(), denX.max(), 10*np.size(domain,1))
        yi = np.linspace(denY.min(), denY.max(), 10*np.size(domain,0))
    # Interpolate for plotting
    zi = griddata((denX, denY), denZ, (xi[None,:], yi[:,None]), method='nearest')
    # range of colorbar by removing data 
    zmin = 0
    zmax = 1
    zi[(zi<zmin) | (zi>zmax)] = None
    # REORGANIZE gridbox INTO HEAT MAP
    domainy = np.repeat(gridbox_y, len(gridbox_x), axis = 0)
    domainx = np.tile(gridbox_x, (len(gridbox_y), 1))
    domainx = np.asarray(list(itertools.chain(*domainx)))
    domainz =  np.asarray(list(itertools.chain(*boxint)))
    # Interpolate for plotting
    boxintplot = griddata((domainx, domainy), domainz, (xi[None,:], yi[:,None]), method='nearest')
    ################################# SMOOTH INTERPOLATION ###################################################
#     # create x-y points to be used in heatmap
#     interpolaxi = np.linspace(denX.min(), denX.max(), np.size(domain,1))
#     interpolayi = np.linspace(denY.min(), denY.max(), np.size(domain,0))
#     # Interpolate for plotting
#     interpolazi = griddata((denX, denY), denZ, (interpolaxi[None,:], interpolayi[:,None]), method='nearest')
#     # range of colorbar by removing data 
#     zmin = 0
#     zmax = 1
#     interpolazi[(interpolazi<zmin) | (interpolazi>zmax)] = None
#     # Interpolate for plotting
#     boxinterpolaplot = griddata((domainx, domainy), domainz, (interpolaxi[None,:], interpolayi[:,None]), method='nearest')
    
    return gridbox_x, gridbox_y, points, xi, yi, boxintplot, boxint, xydframe, raftids, intermedids, nonraftids



##########################################################################################################
##########################################################################################################
############################# UPPER/LOWER LEAFLET DOMAIN AREAS ###########################################
##########################################################################################################
##########################################################################################################
def leafletdomainareas(gridbox_x, gridbox_y, boxint):
    from functions import fillrin,fill2rin
    import copy
    import numpy as np
    # DOMAIN TYPE:
    # 2 = RAFT
    # 1 = INTERMEDIATE
    # 0 = NONRAFT
    for domaintype in [2,1,0]:
        ##########################################################################################
        # LEAFLET DOMAINS' AREAS
        ##########################################################################################
        domainx = gridbox_x
        domainy = gridbox_y 
        domaindx = round(domainx[1]-domainx[0],2)
        domaindy = round(domainy[1]-domainy[0],2)
        dA = domaindy*domaindx

        advnum = 3
        domainarea = []
        #### Lower Raft ####
        boxgridraft = copy.deepcopy(boxint)
        for iy, dy in enumerate(gridbox_y):
            for ix, dx in enumerate(gridbox_x):
                if domaintype == 2:
                    if boxgridraft[iy][ix] == 1:
                        boxgridraft[iy][ix] = 0
                    if boxgridraft[iy][ix] == 2:
                        areaa = 0
                        boxgridraft, area = fillrin(boxgridraft,domaintype,ix,iy,advnum,areaa)
                        domainarea.append(area)
                        advnum += 1
                elif domaintype == 1:
                    if boxgridraft[iy][ix] == 2:
                        boxgridraft[iy][ix] = 0
                    if boxgridraft[iy][ix] == 1:
                        areaa = 0
                        boxgridraft, area = fill2rin(boxgridraft,domaintype,ix,iy,advnum,areaa)
                        domainarea.append(area)
                        advnum += 1
                elif domaintype == 0:
                    if boxgridraft[iy][ix] == 0:
                        areaa = 0
                        boxgridraft, area = fillrin(boxgridraft,domaintype,ix,iy,advnum,areaa)
                        domainarea.append(area)
                        advnum += 1
        if domaintype == 0:
            boxgridraft = [np.where(x==1,0,x) for x in boxgridraft]
            boxgridraft = [np.where(x==2,0,x) for x in boxgridraft]
        # replace (3+) to (3+)-2 = (1+)
        boxgridraft = [np.where(x>=3,x-2,x) for x in boxgridraft]
        # replace 0 to nan
        boxgridraft = [np.where(x==0,np.nan,x) for x in boxgridraft]
        domainarea = list(map(lambda x: x * dA, domainarea))
        totalarea = round((domainx[-1]+domainx[0])*(domainy[-1]+domainy[0]),4)
        domainfromzero = list(map(lambda x: x, boxgridraft))

        if domaintype == 2:
            raftA = np.sum(domainarea)
            raftAcount = np.count_nonzero(domainarea)
            raftdomain = domainarea.copy()
        elif domaintype == 1:
            interA = np.sum(domainarea)
            interAcount = np.count_nonzero(domainarea)
            interdomain = domainarea.copy()
        elif domaintype == 0:
            nonraftA = np.sum(domainarea)
            nonraftAcount = np.count_nonzero(domainarea)
            nonraftdomain = domainarea.copy()

        # REORGANIZE gridbox INTO HEAT MAP FOR CLUSTERS
        # other required variables in function: denX, denY, domain
        
        # number of nearest neighbor interpolation points
        # to expand the size of the domains for visibility
#         intrpltnum = 10
#         xi = np.linspace(denX.min(), denX.max(), intrpltnum*len(domain))
#         yi = np.linspace(denY.min(), denY.max(), intrpltnum*len(domain))
#         domainy = np.repeat(gridbox_y, len(gridbox_x), axis = 0)
#         domainx = np.tile(gridbox_x, (len(gridbox_y), 1))
#         domainx = np.asarray(list(itertools.chain(*domainx)))
#         domainz =  np.asarray(list(itertools.chain(*domainfromzero)))
#         # Interpolate for plotting
#         domainfromzeroplot = griddata((domainx, domainy), domainz, (xi[None,:], yi[:,None]), method='nearest')
    return raftdomain,interdomain,nonraftdomain,raftA,interA,nonraftA,raftAcount,interAcount,nonraftAcount