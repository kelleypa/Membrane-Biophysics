#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2021.31.05 V2 - Patrick Kelley

##########################################################################################################
###################### IMPORT FUNCTIONS ##################################################################
##########################################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import unique

from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.cm as cm
import matplotlib.patches as mpatches

from shapely.geometry import Point, MultiPoint, Polygon

import csv
import re
import time
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.patches import Polygon
import matplotlib.colors as mcolors

#from os import remove, system, path
#from sys import argv, stdout
import subprocess
from sys import argv, stdout
import os
import math

from os import listdir
from os.path import isfile, join
##########################################################################################################
###################### DEFINE SOME FUNCTIONS #############################################################
##########################################################################################################
def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.
    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.
    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

def voronoi_ninebyninebox(points, x_box, y_box):
    point00 = points.copy()
    point10 = points.copy()
    point20 = points.copy()
    point01 = points.copy()
    point02 = points.copy()
    point11 = points.copy()
    point12 = points.copy()
    point21 = points.copy()
    point22 = points.copy()

    point00[:,0] = point00[:,0] - x_box
    point00[:,1] = point00[:,1] - y_box
    point10[:,0] = point10[:,0] - x_box
    point20[:,0] = point20[:,0] - x_box
    point20[:,1] = point20[:,1] + y_box
    point01[:,1] = point01[:,1] - y_box
    point02[:,0] = point02[:,0] + x_box
    point02[:,1] = point02[:,1] - y_box
    point12[:,0] = point12[:,0] + x_box
    point21[:,1] = point21[:,1] + y_box
    point22[:,0] = point22[:,0] + x_box
    point22[:,1] = point22[:,1] + y_box
    
    ninebyninebox = np.concatenate((point00,point10,point20,point01,point02,point11,point12,point21,point22),axis=0)
    
    return ninebyninebox

def voronoi_polygons(points, x_box, y_box):
   
    random_seeds = voronoi_ninebyninebox(points, x_box, y_box)

    vor = Voronoi(random_seeds)
    regions, vertices = voronoi_finite_polygons_2d(vor)

    polygons = []
    # for the middle voronoi tesselled region
    for reg in regions[len(points)*5:len(points)*6]:
        polygon = vertices[reg]
        polygons.append(polygon)
    return polygons

def plot_polygons(polygons, ax=None, alpha=1, linewidth=0.7, saveas=None, show=True):
    from matplotlib.patches import Polygon
    # Configure plot 
    if ax is None:
        plt.figure(figsize=(20,20))
        ax = plt.subplot(111)
    # Remove ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("equal")
    # Set limits
    ax.set_xlim(-x_box,2*x_box)
    ax.set_ylim(-y_box,2*y_box)

    # Add polygons 
    for p, poly in enumerate(polygons):
        colored_cell = Polygon(poly, linewidth=linewidth, alpha=alpha, facecolor=speed[p], edgecolor="black")
        ax.add_patch(colored_cell)

    if not saveas is None:
        plt.savefig(saveas)
    if show:
        plt.show()

    return ax

def listToString(s): 
    # initialize an empty string
    str1 = ""     
    # traverse in the string  
    for ele in s: 
        str1 += ele  
    # return string  
    return str1 

# FOR ORDER PARAMETER SUMMATION PART 
def find_between( s, first, last, lf = 0 ):
    try:
        start = s.index( first ) + len( first ) - lf
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
    
def find_csv_filenames( path_to_dir, suffix=".csv" ):
    filenames = listdir(path_to_dir)
    return [ filename for filename in filenames if filename.endswith( suffix ) ]

# @author: Sajeewa Walimuni Dewage (a.k.a.   Sajeewa Pemasinghe :) )
# https://sajeewasp.com/iteratively-calculate-center-of-mass-gromacs/

# parse a .gro file
# return a list of coordinates for CoM calculation
def beads_read_gro(file, molecule, atoms):
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

# parse a .gro file
# return a list of coordinates for cos^2(theta) calculation
def read_gro(file, atoms):
    line_counter = 0
    number_of_particles = 0
    first, second = [], []
    for line in open(file):
        if line_counter == 1:
            number_of_particles = int(line)
        elif line_counter > 1 and line_counter < number_of_particles + 2:
            if line[10:15].strip() == atoms[0]:
                first.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
            elif line[10:15].strip() == atoms[1]:
                second.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
        line_counter += 1
    return [first, second]
##########################################################################################################
###################### REAL PROGRAM ######################################################################
##########################################################################################################

if len(argv) != 14:
    # coments/usage
    print('''
    Compute the raft vs nonraft phase, defined as:

    Raft = sphingomyelin (DPSM) or cholesterol (CHOL) surrounded within 1nm radius by nothing but DPSM and CHOL lipid centers, with neighboring molecules sharing equal or more than 2/3 (0.66) edges with DPSM/CHOL
    
    NonRaft = anything else

    Currently, the lipids DPSM and CHOL can be analyzed as raft-like with this script.
    It is quite simple to change the definition by editing this program. 
    
    Usage: %s <step1-CoM file> <edr file> <initial time [ns]> <final time [ns]> <skip frames> <dt [picoseconds]> <#lipids> <lipid type>

    > python3 %s my.xtc my.tpr 0 0 1 my.edr  0 6000000 2499 0.03 7000 225 PDPC

    will for example read the frames between 0 and 6,000,000 ps (i.e. 6 microseconds) of the center of mass positions obtained from [VMD -> Extensions -> TK Console] running the step1-CG-CoM.tcl script that outputs my.dat file, md simulation running with time steps of 0.03 picoseconds (30 femtoseconds) saved every 7000 steps, so saved every 0.03*7000 = 210ps, the energy file my.edr to generate PDPC_box_size_0-6mus.xvg file containing the (time, x, y, z, volume), expecting 225 PDPC lipids, determining raft vs nonraft phase for every 2499th frame of each saved frame, so computing every 2499*210ps = 524,790 ps (i.e. 0.525 microseconds). 

    The output is written to folder called:
        PDPC_Phase_Determination/ 
    with files by simulation time {simtime} called: 
        NearestNeighbor_{simtime}mus.csv, Raft_{simtime}mus.csv,and NonRaft_{simtime}mus.csv,
    and with images by simulation time {simtime} called: 
        Voronoi_{simtime}mus.png, Scatter_{simtime}mus.png, and RaftvsNonRaft_{simtime}mus.png
        
    Number of argument inputs is: %i

    ''' % (argv[0], argv[0], len(argv)))
    exit(0) #unsuccessful exit

else:
    start_frame = time.time()
    
    # INPUT
    # snapshots
    trajfile = argv[1]
    tprfile = argv[2]
    # (normalized) orientation of bilayer normal
    orientation_of_bilayer_normal = [float(argv[3]), float(argv[4]), float(argv[5])]
    norm = math.sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
    for i in range(3):
        orientation_of_bilayer_normal[i] /= norm
        stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
        orientation_of_bilayer_normal[0], \
        orientation_of_bilayer_normal[1], \
        orientation_of_bilayer_normal[2]  \
    ))
    stdout.write("\n")
    edrfile = argv[6]
    initial_time0 = int(argv[7])
    final_time0 = int(argv[8])
    vmdstep = int(argv[9])
    dt = float(argv[10])
    nstout = int(argv[11])
    # number of lipids
    number_of_lipids = int(argv[12])
    # lipid type
    phospholipid = argv[13]# NEED ALGORITHM TO READ SYSTEM MOLECULES !!!!
    lipid1 = 'DPSM'# NEED ALGORITHM TO READ SYSTEM MOLECULES !!!!
    lipid2 = 'CHOL'# NEED ALGORITHM TO READ SYSTEM MOLECULES !!!!
    
    
    
    # renaming variables for Step 3
    if phospholipid == 'PDPC':
        phospholipid = 'PUPC'
    traj_skip = vmdstep
    lipid_type = phospholipid
    
    ##########################################################################################################
    ############################## TIME IN FILE CORRECTION  ##################################################
    ##########################################################################################################   
#     initial_time = int(math.ceil(initial_time/(dt*nstout*vmdstep))*(dt*nstout*vmdstep))
#     final_time = int(math.floor(final_time/(dt*nstout*vmdstep))*(dt*nstout*vmdstep))
    initial_time = int(math.ceil(initial_time0/(dt*nstout))*(dt*nstout))
    final_time = range(int(math.ceil(initial_time0/(dt*nstout))*(dt*nstout)),final_time0,int((vmdstep)*(dt)*(nstout)))[-1]
    ########################################################################################################## 
    
    stdout.write("\n initial_time: "+ str(initial_time) + ", " + "final_time: " + str(final_time) + ", " + "traj_skip: "+ str(vmdstep) + "\n\n" )
    
    # Create folder for all analysis
    pathname = './'+phospholipid+'_Analysis_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps/'
    # Check whether the specified pathname exists or not
    isExist = os.path.exists(pathname)
    
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(pathname)
        stdout.write("New analysis directory is created! \n")
    # Checking if the list is empty 
    elif isExist:
        dir = os.listdir(pathname)
        if len(dir) == 0:
            stdout.write("New analysis directory was previously created but empty! Continuing... \n")
        elif len(dir) != 0:
            stdout.write("Analysis folder is not empty. Should we discontinue? Yes:anykey , No:n \n")
            choice = input("> ")    
            if choice == "n":
                stdout.write("Continuing... \n")
            else: 
                stdout.write("Looks like the calculations were already performed. No need to redo it! Closing now...bye \n")
                exit(0) #repetitive exit

    # Create folder for structure files
    structpathname = './'+phospholipid+'_Structure_Files_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps/'
    # Check whether the specified pathname exists or not
    isExist = os.path.exists(structpathname)
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(structpathname)
        stdout.write("New structure (.gro) directory is created! \n\n")

    # Output all frame using trjconv 
    dir = os.listdir(structpathname)
    if len(dir) == 0:
        stdout.write("Output all coordinate files \n")
        #command = "gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %s/frame_dump_.gro > /dev/null" % (trajfile, tprfile, initial_time, final_time, traj_skip,structpathname)
        command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %s/frame_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time, final_time, traj_skip,structpathname)
        print(command)
        subprocess.call(command, shell=True)
    else:
        stdout.write("\nStructure folder is not empty. Should we run gromacs 'gmx trjconv' to find .gro structure files? Yes:anykey , No:n \n")
        choice = input("> ")    
        if choice == "n":
            stdout.write("Continuing... \n\n")
        else: 
            stdout.write("Output all coordinate files \n")
            
            command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o %s/frame_dump_.gro > /dev/null" % ("System", trajfile, tprfile, initial_time, final_time, traj_skip,structpathname)
            print(command)
            subprocess.call(command, shell=True)


    # Output all time, x, y, z, volume of box size from energy .edr file to xgrace .xvg file
    xvgfile = pathname+phospholipid+'_box_size_'+str(initial_time)+'-'+str(final_time)+'ps.xvg'
    
    dir = os.listdir(pathname)
    if len(dir) == 0:
        stdout.write("Output the time, x, y, z, and volume of box size \n")
        command = "gmx energy -f %s -o %s > /dev/null" % (edrfile, xvgfile)
        #command = "echo Box-X Box-Y Box-Z Volume | gmx energy -f %s -o %s > /dev/null" % (edrfile, xvgfile)
        #command = "echo 13 14 15 16 | gmx energy -f %s -o %s > /dev/null" % (edrfile, xvgfile)
        print(command)
        subprocess.call(command, shell=True)
        t, x, y, z, V = [], [], [], [], []
        with open(xvgfile) as f:
            for line in f:
                cols = line.split()
                if cols[0] == '@' or cols[0] == '#':
                    next
                elif len(cols) == 5:
                    t.append(float(cols[0]))
                    x.append(float(cols[1]))
                    y.append(float(cols[2]))
                    z.append(float(cols[3]))
                    V.append(float(cols[4]))
        src = xvgfile
        dst = pathname+phospholipid+'_box_size_'+str(int(t[0]))+'-'+str(int(t[-1]))+'ps.xvg'
        os.rename(src, dst)
    elif len(dir) != 0:
        for fil in dir:
            # search given pattern in the line 
            match = re.search("\.xvg$", fil)
            # if match is found
            if match:
                stdout.write("The file ending with .xvg is: " + fil+'\n')
                xvgmatchedfile = fil
        stdout.write("\nShould we run gromacs 'gmx energy' to find box size and store in xgrace .xvg file anyway? Yes:anykey , No:n \n")
        choice = input("> ")    
        if choice == "n":
            stdout.write("Continuing... \n\n")
            xvgfile = pathname+xvgmatchedfile
            t, x, y, z, V = [], [], [], [], []
            with open(xvgfile) as f:
                for line in f:
                    cols = line.split()
                    if cols[0] == '@' or cols[0] == '#':
                        next
                    elif len(cols) == 5:
                        t.append(float(cols[0]))
                        x.append(float(cols[1]))
                        y.append(float(cols[2]))
                        z.append(float(cols[3]))
                        V.append(float(cols[4]))
        else: 
            stdout.write("Output the time, x, y, z, and volume of box size \n")
            command = "gmx energy -f %s -o %s > /dev/null" % (edrfile, xvgfile)
            #command = "echo 13 14 15 16 | gmx energy -f %s -o %s > /dev/null" % (edrfile, xvgfile)
            print(command)
            subprocess.call(command, shell=True)
            t, x, y, z, V = [], [], [], [], []
            with open(xvgfile) as f:
                for line in f:
                    cols = line.split()
                    if cols[0] == '@' or cols[0] == '#':
                        next
                    elif len(cols) == 5:
                        t.append(float(cols[0]))
                        x.append(float(cols[1]))
                        y.append(float(cols[2]))
                        z.append(float(cols[3]))
                        V.append(float(cols[4]))
            src = xvgfile
            dst = pathname+phospholipid+'_box_size_'+str(int(t[0]))+'-'+str(int(t[-1]))+'ps.xvg'
            os.rename(src, dst)
    
    ##########################################################################################################
    ################################### TOP-BOTTOM LEAFLETS  #################################################
    ##########################################################################################################
    
    leaflet = ["upper","lower"]     
    
    ##########################################################################################################
    ############################## CALCULATE CoM and [cos(theta)]^2  #########################################
    ##########################################################################################################   

    comstrFile = pathname+phospholipid+'_CoM_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps.dat'
    isExist = os.path.exists(comstrFile)
    
    # FUTURE: modify this to be automatically generated as input argv[] 
    LipidNames = [phospholipid, "DPSM", "CHOL"] # NEED ALGORITHM TO READ SYSTEM MOLECULES !!!!
    
    ########################### BEADS, BACKBONEBEADS, AND MASS #############################################

    if phospholipid == 'POPC':
        Beads = [["C1A", "D2A", "C3A", "C4A", "GL1", "GL2", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"],
                ["AM2", "AM1", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"], 
                ["C2", "C1", "R5", "R4", "R3", "R2", "R1", "ROH"]]
        Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72],
                [72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72], 
                [72, 39.49, 0, 0, 159.65, 38.69, 0, 77.22]]           
        BackBoneBeads = [["GL1", "GL2"], 
                        ["AM2", "AM1"],
                        ["ROH"]]
    elif phospholipid == 'PUPC':
        Beads = [["D1A", "D2A", "D3A", "D4A", "D5A", "GL1", "GL2", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"],
                ["AM2", "AM1", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"], 
                ["C2", "C1", "R5", "R4", "R3", "R2", "R1", "ROH"]]        
        Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72],
                [72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72], 
                [72, 39.49, 0, 0, 159.65, 38.69, 0, 77.22]]
        BackBoneBeads = [["GL1", "GL2"], 
                        ["AM2", "AM1"],
                        ["ROH"]]
    if phospholipid == 'POPE':
        Beads = [["C1A", "C2A", "C3A", "C4A", "GL1", "GL2", "C1B", "C2B", "D3B", "C4B", "C5B", "NH3", "PO4"],
                ["AM2", "AM1", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"], 
                ["C2", "C1", "R5", "R4", "R3", "R2", "R1", "ROH"]]
        Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72],
                [72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72], 
                [72, 39.49, 0, 0, 159.65, 38.69, 0, 77.22]]           
        BackBoneBeads = [["GL1", "GL2"], 
                        ["AM2", "AM1"],
                        ["ROH"]]
    elif phospholipid == 'PUPE':
        Beads = [["D1A", "D2A", "D3A", "D4A", "D5A", "GL1", "GL2", "C1B", "C2B", "C3B", "C4B", "NH3", "PO4"],
                ["AM2", "AM1", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"], 
                ["C2", "C1", "R5", "R4", "R3", "R2", "R1", "ROH"]]        
        Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72],
                [72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72], 
                [72, 39.49, 0, 0, 159.65, 38.69, 0, 77.22]]
        BackBoneBeads = [["GL1", "GL2"], 
                        ["AM2", "AM1"],
                        ["ROH"]]
    # Create Dictionary with Dataframes of Beads' Mass (easier for BBB lookup in calculation)
    df_beadmass = {}
    for l,lip in enumerate(LipidNames):
        d = {}
        d = {Beads[l][b]: [Mass[l][b]] for b in range(len(Beads[l]))}
        df_beadmass[lip] = pd.DataFrame(data=d)
    ######################################## BONDS ##########################################################
    
    phosphatidylcholine_bond_names = "NC3-GL1 "
    phosphatidylethanolamine_bond_names = "NH3-GL1 "
    # PCs
    if   lipid_type == "POPC": bond_names = phosphatidylcholine_bond_names + "PO4-C1A PO4-GL2 GL1-D2A C1A-C3A D2A-C4A GL2-C2B C1B-C3B C2B-C4B\n"
    elif lipid_type == "PUPC": bond_names = phosphatidylcholine_bond_names + "PO4-D1A PO4-GL2 GL1-D2A D1A-D3A D2A-D4A D3A-D5A GL2-C2B C1B-C3B C2B-C4B\n"
        
    # PEs
    elif lipid_type == "PUPE": bond_names = phosphatidylethanolamine_bond_names + "PO4-D1A PO4-GL2 GL1-D2A D1A-D3A D2A-D4A D3A-D5A GL2-C2B C1B-C3B C2B-C4B\n"
    elif lipid_type == "POPE": bond_names = phosphatidylethanolamine_bond_names + "PO4-C1A PO4-GL2 GL1-C2A C1A-C3A C2A-C4A GL2-C2B C1B-D3B C2B-C4B D3B-C5B\n"
    
    # DPSM
    dpsm_bond_names = "NC3-AM1 AM1-C2A AM2-C2B T1A-C3A C1B-C3B C2B-C4B\n"
    # CHOL (newer V2_02 itp version)
    chol_bond_names = "ROH-C1 R3-C2 R2-C2\n"
    
    bond_list = [bond_names, dpsm_bond_names, chol_bond_names]
    ##########################################################################################################
        

    if BackBoneBeads != Beads:
        comstrFile = pathname+phospholipid+'_BBB_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps.dat'
        isExist = os.path.exists(comstrFile)

    if not isExist:
        # Nested Dictionary of Bonds dictionaries inside Lipids dictionaries
        df_bond = {}
        ###################### center of mass PARAMETERs #########################################################
        com_count = 0
        df_com = pd.DataFrame(data={"time": [], "lipid": [], "id": [], "resid_start": [], "resid_end": [], "leaflet": [],  "x": [], "y": [], "z": []})
        # convert column "time", "id", "resid_start", and "resid_end" to int64 dtype
        df_com = df_com.astype({"time": int, "id": int, "resid_start": int, "resid_end": int})
        ##########################################################################################################
        stdout.write("Starting center of mass calculation & [cos(theta)]^2 calculation \n\n")
        ##########################################################################################################
        ################## For each dumped frame (aka timesnap) ###################################################
        ##########################################################################################################
        file_count = 0
        while True:
            filename = structpathname+"frame_dump_" + str(file_count) + ".gro"
            if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
                break

            tempora = initial_time + (file_count)*(vmdstep)*(dt)*(nstout)
            stdout.write("\n Taking care of snapshot %s for time %i ps \n" % (filename,tempora))
            
            # Find all bead locations for all lipids and store dataframes of positions in dictionary
            # Python Dictionary titled by LipidNames
            df_file = {}
            for n, name in enumerate(LipidNames):
                df_xyz = beads_read_gro(filename,name,Beads[n])
                df_file[name] = df_xyz
            
            ############## CoM for each lipid in BackBoneBeads list CALCULATION ######################################
            stdout.write("Starting center of mass calculation for all specified beads in each lipid \n")
            
            for mol1,molecul1 in enumerate(LipidNames):
                df_xyz = df_file[molecul1]
                for lip in range(len(df_xyz)):
                    # vector sum of all beads designated for each lipid (row of dataframe)
                    vector = [0.0, 0.0, 0.0]
                    mass_sum = 0
                    for bbb in BackBoneBeads[mol1]:
                        vector[0] = vector[0] + float(df_beadmass[molecul1][bbb])*df_xyz.loc[lip, 'x_'+bbb]
                        vector[1] = vector[1] + float(df_beadmass[molecul1][bbb])*df_xyz.loc[lip, 'y_'+bbb]
                        vector[2] = vector[2] + float(df_beadmass[molecul1][bbb])*df_xyz.loc[lip, 'z_'+bbb]
                        mass_sum = mass_sum + float(df_beadmass[molecul1][bbb])
                    # center of mass of selected beads [assume beads are all weight: 72 amu...except CHOL]
                    # source for assumption: http://www.ks.uiuc.edu/Training/Tutorials/martini/rbcg-tutorial.pdf pg 11.
                    comvector = [x/mass_sum for x in vector]
                    df_com.loc[com_count, ['x','y','z']] = comvector
                    df_com.loc[com_count, ['time','lipid','id','resid_start', 'resid_end']] = [int(tempora), molecul1, int(df_xyz.loc[lip,'id']), int(df_xyz.loc[lip,'resid_start']), int(df_xyz.loc[lip,'resid_end'])]
                    com_count += 1
            ######################### Upper and Lower Leaflet Determination ########################################
            avg_z = df_com['z'].loc[df_com['time'] == tempora].mean()
            for lip in df_com['z'].loc[df_com['time'] == tempora].index:
                if df_com.loc[lip,'z'] >= avg_z:
                    df_com.loc[lip, 'leaflet'] = 'upper'
                elif df_com.loc[lip,'z'] < avg_z:
                    df_com.loc[lip,'leaflet'] = 'lower'
            ###################### (cos(theta))^2 PARAMETERs #########################################################
            if file_count == 0:
                for mol,molecul in enumerate(LipidNames):
                    xydframe = df_com.loc[df_com['lipid'] == molecul] # need df_com info to build df_bond
                    df_bond[molecul] = {}
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
            ###################### (cos(theta))^2 CALCULATION #########################################################
            for mol,molecul in enumerate(LipidNames):
                stdout.write("Starting [cos(theta)]^2 calculation for all specified bonds in %s \n" % (molecul))
                # compute order parameter for each bond, for each snapshot
                df_xyz = df_file[molecul]
                bonds = bond_list[mol].split()
                for bond in df_bond[molecul]:
                    df_bond[molecul][bond].loc[file_count,'time'] = tempora
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
                        norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
                        # compute projection on the bilayer normal
                        projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
                        # order parameter
                        order_parameter = projection**2/norm2# convert strings to floats
                        # Append to dataframe in dictionary
                        df_bond[molecul][bon].loc[file_count,str(lipid_id)] = order_parameter
            file_count += 1
        
        # write the center of mass results to csv file
        stdout.write("\n The center of mass of the selected lipids has been computed. Saving data to analysis folder. \n\n")
        df_com.to_csv(comstrFile, sep=',',index=False)
        
        # write the [cos(theta)]^2 for each bond to csv file
        for mol, molecul in enumerate(LipidNames):
            bonds = bond_list[mol].split()
            for bond in bonds:
                datafilename = pathname+molecul+'_cos2_'+bond+'.dat'
                df_bond[molecul][bond].to_csv(datafilename, sep=',', index=False)
    
    else:
        stdout.write("The center of mass file already exists. Importing data. \n\n")
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
    

    ##########################################################################################################
    ###################### START PHASE DETERMINATION #########################################################
    ##########################################################################################################
    stdout.write("Do you want plots with raft, INTERMEDIATE and nonraft? Yes:y , No:anykey \n")
    choice = input("> ")    
    if choice == "y":
        stdout.write("Will plot raft, intermediate, and nonraft. \n")
        intermed = True
    else: 
        stdout.write("Will plot raft and nonraft only. \n\n")
        intermed = False
    
    
        
    
   
    
    # MAIN LOOP
#     for fram, tempora in enumerate(range(int(math.ceil(initial_time/(dt*nstout))*(dt*nstout)),final_time,int((vmdstep)*(dt)*(nstout)))):
    
    for leaf in leaflet:
        # INITIALIZING GLOBAL RAFT & NONRAFT INDICE LISTS
        raftindexfile = pathname+leaf+'_'+phospholipid+'_raft-nonraft-indices_'+str(initial_time)+'-'+str(int(vmdstep*dt*nstout))+'-'+str(final_time)+'ps'
        raftlist = [[]]*len(df_com.time.unique())
        nonraftlist = [[]]*len(df_com.time.unique())
        timelist = [[]]*len(df_com.time.unique())
        
        # Anounce which leaflet we are working on
        stdout.write("\n%s\n   Starting %s bilayer leaflet. \n%s\n" % ("*"*46,leaf,"*"*46))
        
        for fram, tempora in enumerate(range(initial_time,final_time0, int((vmdstep)*(dt)*(nstout)))):
            frame = str(tempora)
            ########## TIME #########################################################################################
            simtime = tempora*1e-6 # microseconds
            timelist[int(fram)] = tempora # picoseconds
            ########################################################################################################
            #df_xy = {} # make dictionary for xy positions of all lipids in upper and lower leaflet
                
            xydframe = df_com.loc[df_com['leaflet'] == leaf] # coordinates of all lipids in select frame
            xydframe = xydframe.loc[xydframe['time'] == tempora]

            xydframe.index = range(len(xydframe)) # set the indices to 0 of upper or lower leaflet 
            
            df_xyd = xydframe.loc[:,['id','lipid', 'resid_start','resid_end']] # contains **id**, lipid, resid_start, resid_end for phase characteristics
            xydframe = xydframe.loc[:,['id','x','y']] # only contains **id** and x, y coordinate for phase determination
        
            ##################### PLOT COLOR OF PHOSPHOLIPID, CHOLESTEROL, and SPHINGOMYELIN #########################
            # PHOSPHOLIPID
            speed_LIPID = ['b' for _ in range(len(df_xyd.loc[df_xyd['lipid']== phospholipid]))]
            # DPSM
            speed_SM = ['r' for _ in range(len(df_xyd.loc[df_xyd['lipid']== 'DPSM']))]
            # CHOL
            speed_CHOL = ['w' for _ in range(len(df_xyd.loc[df_xyd['lipid']== 'CHOL']))]
            # Color scheme for all molecules
            speed = np.concatenate((speed_LIPID, speed_SM, speed_CHOL))
            speed = np.concatenate((speed,speed,speed,speed,speed,speed,speed,speed,speed),axis=0)
            ##########################################################################################################            
            # Find box size from imported values of the Energy File my.edr
            timeindex = t.index(tempora)
            x_box = x[timeindex] #nm 
            y_box = y[timeindex] #nm 
            
            stdout.write('\n%i ps, X_Box = %.2f nm, Y_Box = %.2f nm \n' % (tempora,x_box,y_box))

            ##########################################################################################################            
            from shapely.geometry import Point, MultiPoint, Polygon
            points = xydframe.loc[:,['x','y']].to_numpy()
            
            # compute Voronoi tesselation
            polygons = voronoi_polygons(points, x_box, y_box)
            #plot_polygons(polygons, ax=None, alpha=1, linewidth=0.7, saveas=None, show=False)
            
            ##########################################################################################################
            ################################ NEAREST NEIGHBOR: Within 1nm Radius #####################################
            ##########################################################################################################
            strFile = pathname+leaf+'_nearestneighbor_'+phospholipid+'_'+frame+'ps.csv'
            # Check whether the specified pathname exists or not
            isExist = os.path.exists(strFile)
            if not isExist:

                # starting time
                start = time.time()

                # *so*rted *po*sition of lipids within Distance of each other
                sopo=xydframe.copy()
                #sopo=sopo.rename(columns={0: "X1", 1:"X2"})
        
                box_list = [[]]*len(sopo)
                nm = 1 # radius to look around for nearest neighbors
                for p in range(0,len(sopo)):
                    for a in range(p+1,len(sopo)):
                        if math.sqrt((sopo.loc[p,'x']-sopo.loc[a,'x'])**2+(sopo.loc[p,'y']-sopo.loc[a,'y'])**2) <= nm:
                            #box_list[p] = box_list[p] + [sopo.iloc[a].name]
                            #box_list[a] = box_list[a] + [sopo.iloc[p].name] #can't use row name as index, due to flip-flops
                            box_list[p] = box_list[p] + [int(sopo.loc[a,'id'])]
                            box_list[a] = box_list[a] + [int(sopo.loc[p,'id'])]
                # Box List holds all nearest neighbor molecules withn nm distance 
                # INDEX 0-224: LIPIDS, INDEX 225-449: DPSM, INDEX 450-674: CHOL
            
                # Create a new file
                df = pd.DataFrame(data={"NN": box_list})
                df.to_csv(strFile, sep=',',index=False)

                stdout.write("Nearest Neighbor List was computed and saved. The new file is created in folder.")

                # end time
                end = time.time()
                # total time taken
                stdout.write("Runtime of the program is %.2f second \n" % (end-start))    
            else:
                stdout.write("The file already exists. Importing data.")
                test_box_list =pd.read_csv(strFile,index_col=False) 
                test_list = test_box_list.values.tolist()
                box_list = []
                for i in test_list:
                    tt = listToString(i)
                    test = re.split('\[|, |\]|\n',tt)
                    
                    test = ' '.join(test).split()
                    test = list(map(int, test))
                    box_list.append(test)
                    
            # nearest neighbor (column: nn) - box id (column: id) dataframe 
            df_box = pd.DataFrame(xydframe['id'])
            df_box['nn'] = box_list
            # highest SM/CHOL density calculation 
            box_cnt = []
            for i, ind in enumerate(box_list):
                cnt = 0
                for ii, nnR in enumerate(ind):
                    lipid_typ = df_xyd.loc[df_xyd['id'] == nnR,'lipid'].item()
                    if lipid_typ == 'CHOL' or lipid_typ == 'DPSM':
                        cnt += 1
                lngth = len(ind)
                if lngth != 0:
                    box_cnt.append([cnt/lngth,cnt,lngth])
                else: box_cnt.append([0,cnt,lngth])
            ##########################################################################################################
            ############ PHASE DETERMINATION: Cluster Center AS LIPID SURROUNDED BY 100% DPSM OR CHOL ################
            ##########################################################################################################
            start = time.time()
            tessels_list = np.zeros(len(box_list))
            tessels_int_list = np.zeros(len(box_list))
            
            # list of nearest neighbor clusters
            nn_cluster_list = [[]]*len(box_list)
            
            # Nearest Neighbor sphingomyelin & cholesterol threshold
            SMCHOLEdgeThres = 1
            # list of *N*earest *N*eighbor *SUR*rounding lipids already *PROC*cessed
            NNsurproc = []
            # list of lipids that we need to check (calculate) around it to see for other raft-like molecules (set by SMCHOLEdgeThres)
            NNcalc = []
            for i in range(0,len(df_box)):
                SMCHOLdensity = box_cnt[i][0]
                
                if SMCHOLdensity == 1 and tessels_list[i] != 1:
                    
                    # That lipid completely surrounded by CHOL/SM in 1nm radius is cluster center
                    tessels_list[i] = 1
                    # Find all lipids [DSPM/CHOL] surrounding that lipid and are now part of cluster
                    for ii, ind in enumerate(box_list[i]):
                        tessels_list[df_box.loc[df_box['id']==ind].index.item()] = 1
                        nn_cluster_list[i] = nn_cluster_list[i]+[ind]

                    NNsurproc.append(i)
                    
                    # start calculating all surrounding lipids if part of cluster
                    NNcalc = nn_cluster_list[i].copy()
                    while len(NNcalc) != 0:
                        
                        ind = NNcalc[0]
                        
                        if ind not in NNsurproc:
                            for k, iind in enumerate(df_box.loc[df_box['id']==ind, 'nn'].item()):
                                nn_vertex = polygons[df_box.loc[df_box['id'] == iind].index.item()]
                                
                                lngth = len(polygons[df_box.loc[df_box['id']==iind].index.item()])
                                # Count all nearest neighbor DPSM/CHOL lipids that share an edge with the select lipid
                                cnt = 0
                                cnt_ind = 0
                                cnt_choldpsm = 0
                                for j,iiind in enumerate(df_box.loc[df_box['id']==iind, 'nn'].item()):
                                    cnt_ind += 1
                                    lipid_typ = df_xyd.loc[df_xyd['id'] == iiind,'lipid'].item()
                                    if lipid_typ == 'CHOL' or lipid_typ == 'DPSM':
                                        cnt_choldpsm +=1
                                        for b, vertab in enumerate(nn_vertex):
                                            for j, vertij in enumerate(polygons[df_box.loc[df_box['id']==iiind].index.item()]):
                                                if np.allclose(vertab,vertij):            
                                                    cnt += 1
                                
                                #print(nn_vertex)
                                if cnt_choldpsm != 0:
                                    rashio = cnt/(lngth*cnt_choldpsm)
                                else: rashio = 0
                                #print(iind, cnt_ind, cnt_choldpsm, cnt, lngth, cnt/(lngth), rashio)
                                # Check if edges above SMCHOLEdgeThres
                                #cnt = cnt/6.6667 # double counting happening
                                if rashio >= SMCHOLEdgeThres:
                                    tessels_list[df_box.loc[df_box['id']==iind].index.item()] = 1
                                    NNcalc.append(iind)
                                    nn_cluster_list[i] = nn_cluster_list[i]+[iind]
                                elif 0.5 <= rashio < 1:
                                    tessels_int_list[df_box.loc[df_box['id']==iind].index.item()] = 1
                                else: NNsurproc.append(iind)  
                            NNsurproc.append(ind)
                        NNcalc.pop(0)

            #Remove clusters less than certain size
            for i in range(len(nn_cluster_list)):
                cnt = np.count_nonzero(nn_cluster_list[i])
                if not cnt == 0 and cnt < 4:
                    tessels_list[i] = 0
                    for ind in nn_cluster_list[i]:
                        tessels_list[df_box.loc[df_box['id']==ind].index.item()] = 0

            end = time.time()
            stdout.write("Runtime of the Phase Determination is %.2f second \n" % (end-start)) 
            ##########################################################################################################
            ######################### SAVE Raft vs NonRaft positions #################################################        
            ##########################################################################################################
            raftstrFile = pathname+leaf+'_raft_'+phospholipid+'_'+frame+'ps.csv'
            nonraftstrFile = pathname+leaf+'_nonraft_'+phospholipid+'_'+frame+'ps.csv'
            # Check whether the specified pathname exists or not
            isExist = os.path.exists(raftstrFile)
            if not isExist:
                raft_indices = []
                nonraft_indices = []

                #     raft_indices = [i for i,t in enumerate(tessels_list) if t==1]
                for i,tt in enumerate(tessels_list):
                    lipid_id = int(df_xyd.loc[i, 'id'])
                    residue_start = int(df_xyd.loc[i, 'resid_start'])
                    residue_end = int(df_xyd.loc[i, 'resid_end'])
                    # Raft-Like
                    if tt==1:
                        raft_indices.append([residue_start, residue_end, lipid_id])
                    #NonRaft-Like
                    else:
                        nonraft_indices.append([residue_start, residue_end, lipid_id])
                raft_dataframe=pd.DataFrame(raft_indices)
                nonraft_dataframe=pd.DataFrame(nonraft_indices)
                # Global Raft and NonRaft List
                raftlist[fram] = raftlist[fram]+raft_indices
                nonraftlist[fram] = nonraftlist[fram]+nonraft_indices
                # Create a new file
                raft_dataframe.to_csv(raftstrFile, header=False, index=False, sep=',', mode='a')
                nonraft_dataframe.to_csv(nonraftstrFile, header=False, index=False, sep=',', mode='a')
                stdout.write("Raft & NonRaft Lipid Indices were saved. The new file is created in folder.")
            else: stdout.write("Files for Raft & NonRaft Lipid Indices already were saved. \n")
            ##########################################################################################################
            ######################### Final Phase Determination Computational Output #################################        
            ##########################################################################################################
            end_frame = time.time()
            stdout.write("Number of raft-like lipids discovered: "+str(np.sum(tessels_list))+"\n")
            stdout.write("Runtime of the entire program so far is %0.2f seconds \n" % (end_frame - start_frame))
            
            
            ##########################################################################################################
            ##########################################################################################################
            ######### SAVE IMAGES: Voronoi Tessellation, Scatter, and Phase Separated Voronoi Tessellation ###########        
            ##########################################################################################################
            ##########################################################################################################
            from matplotlib.patches import Polygon
            # SOME COLORS
            purple = "#7E1E9C"
            violet = "#9A0EEA"
            darkblue = "#030764"
            # Black box around legend colors
            # SOURCE:
            # https://matplotlib.org/3.1.1/tutorials/intermediate/legend_guide.html
            class AnyObject(object):
                pass
            class AnyObjectHandler(object):
                def legend_artist(self, legend, orig_handle, fontsize, handlebox,):
                    global c
                    c += 1
                    colorwork = ['r','w','b','yellow',darkblue]
                    x0, y0 = handlebox.xdescent, handlebox.ydescent
                    width, height = handlebox.width, handlebox.height
                    patch = mpatches.Rectangle([x0, y0], width, height, facecolor=colorwork[c],
                                            edgecolor='black', linewidth=3,
                                            transform=handlebox.get_transform())
                    handlebox.add_artist(patch)
                    return patch
            class AnyObject2(object):
                pass
            class AnyObjectHandler2(object):
                def legend_artist(self, legend, orig_handle, fontsize, handlebox,):
                    global c
                    c += 1
                    colorwork = ['yellow','gray']
                    x0, y0 = handlebox.xdescent, handlebox.ydescent
                    width, height = handlebox.width, handlebox.height
                    patch = mpatches.Rectangle([x0, y0], width, height, facecolor=colorwork[c],
                                            edgecolor='black', linewidth=3,
                                            transform=handlebox.get_transform())
                    handlebox.add_artist(patch)
                    return patch 
            class AnyObject3(object):
                pass
            class AnyObjectHandler3(object):
                def legend_artist(self, legend, orig_handle, fontsize, handlebox,):
                    global c
                    c += 1
                    colorwork = ['yellow','green','gray']
                    x0, y0 = handlebox.xdescent, handlebox.ydescent
                    width, height = handlebox.width, handlebox.height
                    patch = mpatches.Rectangle([x0, y0], width, height, facecolor=colorwork[c],
                                            edgecolor='black', linewidth=3,
                                            transform=handlebox.get_transform())
                    handlebox.add_artist(patch)
                    return patch
            scalebar = ScaleBar(
                1e-9,
                "m",
                length_fraction=0.10,
                rotation="horizontal",
                scale_loc="right",
                border_pad=0.5,
                pad=0.1,
            )
            ##########################################################################################################
            ################### SAVE Voronoi Tessellation Image ######################################################
            ##########################################################################################################
            # Resize plot
            plt.rcParams["figure.figsize"] = (10,10)
            
            # Configure plot 
            ax = plt.subplot(111)

            # Remove ticks
            ax.set_xticks([])
            ax.set_yticks([])

            ax.axis("equal")

            alpha=1
            linewidth=4

            # Set limits
            ax.set_xlim(-2,x_box+2)
            ax.set_ylim(-2,y_box+2)
            ax.add_artist(scalebar)

            rangefrom = 0
            rangeto = len(df_xyd)
            ranger = range(rangefrom,rangeto)
            for rr, r in enumerate(ranger):
                # Add polygons 
                colored_cell = Polygon(polygons[r][:],
                                    linewidth=linewidth, 
                                    alpha=alpha,
                                    facecolor=speed[r],
                                    edgecolor="black")
                ax.add_patch(colored_cell)
            # Creating legend with black box around colors
            c = -1
            pop_a = 'DPSM'
            pop_b = 'CHOL'
            pop_c = phospholipid
            plt.legend([AnyObject(), AnyObject(), AnyObject()], [pop_a, pop_b, pop_c],
                    handler_map={AnyObject: AnyObjectHandler()},prop={'size': 20},bbox_to_anchor=(0,-0.08,1,0.2), loc="lower left",
                            mode="expand", borderaxespad=0, ncol=3)
            plt.title('Lipid: '+phospholipid+'; t='+str(tempora)+'ps',fontsize=35)

            pltfilename = pathname+leaf+'_voronoitessellation_'+phospholipid+'_'+frame+'ps.png'
            plt.savefig(pltfilename)
            plt.cla()
            ##########################################################################################################
            ################## SAVE Phase Separated Voronoi Tessellation Image #######################################
            ##########################################################################################################
        
            ax = plt.gca()
            # Set limits
            ax.set_xlim(-2,x_box+2)

            ax.set_ylim(-2,y_box+2)
            #hide x-axis
            ax.get_xaxis().set_visible(False)
            #hide y-axis
            ax.get_yaxis().set_visible(False)
            ax.add_artist(scalebar)

            plt.scatter(points[:, 0], 
            points[:, 1],c = speed[0:len(df_xyd)],edgecolors='black', linewidth=3, s=200)
            c = -1
            pop_a = 'DPSM'
            pop_b = 'CHOL'
            pop_c = phospholipid
            # plt.legend(handles=[pop_a,pop_b,pop_c,pop_d,pop_e],prop={'size': 20},bbox_to_anchor=(0,-0.08,1,0.2), loc="lower left",
            #                 mode="expand", borderaxespad=0, ncol=3)
            plt.legend([AnyObject(), AnyObject(), AnyObject()], [pop_a, pop_b, pop_c],
                    handler_map={AnyObject: AnyObjectHandler()},prop={'size': 20},bbox_to_anchor=(0,-0.08,1,0.2), loc="lower left",
                            mode="expand", borderaxespad=0, ncol=3)
            plt.title('Lipid: '+phospholipid+'; t='+str(tempora)+'ps',fontsize=35)
            pltfilename = pathname+leaf+'_scatterplot_'+phospholipid+'_'+frame+'ps.png'
            plt.savefig(pltfilename)
            plt.cla()
            ##########################################################################################################
            ############################ SAVE ScatterPlot Image ######################################################
            ##########################################################################################################
            # Threshold values, i.e. how many molecules it needs to have contact with
            nnSMCHOL = 0.66
            nnLIPID = 0.8
            
            # Resize plot
            plt.rcParams["figure.figsize"] = (10,10)

            rangefrom = 0
            rangeto = len(df_xyd)
            ranger = range(rangefrom,rangeto)
            for rr, r in enumerate(ranger):
                poly = polygons[r][:]
                plt.fill(*zip(*poly), color="gray", edgecolor = 'black',linewidth=4)
                
            ax = plt.gca()
            # Set limits
            ax.set_xlim(-2,x_box+2)
            ax.set_ylim(-2,y_box+2)
            #hide x-axis
            ax.get_xaxis().set_visible(False)
            #hide y-axis
            ax.get_yaxis().set_visible(False)
            ax.add_artist(scalebar)

            # Creating legend with black box around colors
            c = -1
            pop_d = 'Raft'
            pop_e = 'Interface'
            pop_f = 'Nonraft'

            

            if intermed == True:
                # INTERFACE REGIONS: Nearest Neighbor-Highest SM/CHOL Density + SM/CHOL Connectivity Threshold    
                rangefrom = 0
                rangeto = len(df_xyd)
                ranger = range(rangefrom,rangeto)
                # plt.plot(points[rangefrom:rangeto,0], points[rangefrom:rangeto,1], 'bo', markersize = 10)
                for rr, r in enumerate(ranger):
                    if tessels_int_list[rr] == 1:
                        poly = polygons[r][:]
                        plt.fill(*zip(*poly), color="green", edgecolor = 'black',linewidth=4)
                plt.legend([AnyObject3(), AnyObject3(), AnyObject3()], [pop_d, pop_e, pop_f],
                        handler_map={AnyObject3: AnyObjectHandler3()},prop={'size': 20},bbox_to_anchor=(0,-0.08,1,0.2), loc="lower left",
                                mode="expand", borderaxespad=0, ncol=3)
            elif intermed == False:
                plt.legend([AnyObject2(), AnyObject2()], [pop_d, pop_f],
                        handler_map={AnyObject2: AnyObjectHandler2()},prop={'size': 20},bbox_to_anchor=(0,-0.08,1,0.2), loc="lower left",
                                mode="expand", borderaxespad=0, ncol=2)
            # RAFT-LIKE REGIONS: Nearest Neighbor-Highest SM/CHOL Density + SM/CHOL Connectivity Threshold    
            rangefrom = 0
            rangeto = len(df_xyd)
            ranger = range(rangefrom,rangeto)
            # plt.plot(points[rangefrom:rangeto,0], points[rangefrom:rangeto,1], 'bo', markersize = 10)
            for rr, r in enumerate(ranger):
                if tessels_list[rr] == 1:
                    poly = polygons[r][:]
                    plt.fill(*zip(*poly), color="yellow", edgecolor = 'black',linewidth=4)  
            #plt.title('Lipid: '+phospholipid+'; Nearest Neighbor SM/CHOL$\geq$'+str(SMCHOLEdgeThres)+'; t='+str(tempora)+'ps',fontsize=20)
            plt.title('Lipid: '+phospholipid+'; t='+str(tempora)+'ps',fontsize=20)
            pltfilename = pathname+leaf+'_raft-nonraftplot_'+phospholipid+'_'+frame+'ps.png'
            plt.savefig(pltfilename)
            plt.cla()
            
        ##########################################################################################################
        ######################### SAVE GLOBAL RAFT & NONRAFT INDICES #############################################        
        ##########################################################################################################
        df_raft = pd.DataFrame(data={'time[ps]': timelist , 'raft': raftlist, 'nonraft': nonraftlist})
        df_raft.to_csv(raftindexfile, sep=',',index=False)
        #stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
        stdout.write("\n Analysis of %s phase domains done.\n " % (leaf))
    ############################### END OF PHASE DETERMINATION ###############################################    
    
    

    
    ##########################################################################################################
    ######################### CALCULATE ORDER PARAMETER  #####################################################        
    ##########################################################################################################
    
    # INDEX REMINDER (DPSM <-> PHOSPHOLIPID: VMD)
    
    # PREVIOUS VERSION:
    # DPSM Indices: 0-224
    # PHOSPHOLIPID Indices: 225-449
    # CHOL Indices: 450-674
    
    # MORE GENERAL VERSION:
    # DSPM Indices:
    # PHOSPHOLIPID Indices: 
    # CHOL Indices: 
    
    for leaf in leaflet: 
        for mol, molecul in enumerate(LipidNames):
            stdout.write("\n Finally summing [cos(theta)]^2 for %s %s - raft and nonraft.\n" % (leaf, molecul))
            filenames = find_csv_filenames(pathname)
            strt = True
            bonds = bond_list[mol].split()
            
            for bond_name in bonds:
                file1 = open(pathname+molecul+'_cos2_'+bond_name+'.dat','r')
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
                # INITIALIZE
                if strt == True:
                    df_raft_op = pd.DataFrame(data=df['time'])
                    df_nonraft_op = pd.DataFrame(data=df['time'])
                    df_raft_op.loc[len(df_raft_op)+1,'time'] = 'Sum'
                    df_nonraft_op.loc[len(df_nonraft_op)+1,'time'] = 'Sum'
                    strt = False

                for indtim, timmie in enumerate(df['time']):
                    for name in filenames:
                        if find_between(name,leaf+"_","_"+phospholipid) == 'raft' and find_between(name,phospholipid+"_","ps.csv") == str(int(timmie)):           
                            raftie = pd.read_csv(pathname+name, header=None)
                            raftiesum = 0
                            raft_cnt = 0
                            for ndx in raftie[2]:
                                if str(ndx) in df_id:
                                    raftiesum = raftiesum + df.loc[indtim,str(ndx)]
                                    raft_cnt += 1
                            df_raft_op.loc[indtim, 'N'] = raft_cnt
                            df_raft_op.loc[indtim, bond_name] = raftiesum

                        elif find_between(name,leaf+"_","_"+phospholipid) == 'nonraft' and find_between(name,phospholipid+"_","ps.csv") == str(int(timmie)):           
                            nonraftie = pd.read_csv(pathname+name, header=None)
                            nonraftiesum = 0
                            nonraft_cnt = 0
                            for ndx in nonraftie[2]:
                                if str(ndx) in df_id:
                                    nonraftiesum = nonraftiesum + df.loc[indtim,str(ndx)]
                                    nonraft_cnt += 1
                            df_nonraft_op.loc[indtim, 'N'] = nonraft_cnt
                            df_nonraft_op.loc[indtim, bond_name] = nonraftiesum

                df_raft_op.loc[df_raft_op['time'] == 'Sum', bond_name] = df_raft_op[bond_name].sum()
                df_nonraft_op.loc[df_nonraft_op['time'] == 'Sum', bond_name] = df_nonraft_op[bond_name].sum()
            df_raft_op.loc[df_raft_op['time'] == 'Sum', 'N'] = df_raft_op['N'].sum()
            df_nonraft_op.loc[df_nonraft_op['time'] == 'Sum', 'N'] = df_nonraft_op['N'].sum()

            # SAVE TO CSV
            opraftpathname = pathname + leaf+'_raft_'+molecul+'_'+str(int(initial_time))+'-'+str(int(final_time))+'ps_OrderParameter.csv'
            df_raft_op.to_csv(opraftpathname, header=True, index=False, sep=',', mode='a')

            opnonraftpathname = pathname + leaf + '_nonraft_'+molecul+'_'+str(int(initial_time))+'-'+str(int(final_time))+'ps_OrderParameter.csv'
            df_nonraft_op.to_csv(opnonraftpathname, header=True, index=False, sep=',', mode='a')

    stdout.write("\n\n Finished, please come again. %s\n" % (" "*56))
