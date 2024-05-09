#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2022.02.28 V7 generalize for all lipids & combine all analysis - Patrick Kelley



def lipidinventory(LipidNames):
    import pandas as pd
    ######################################## BONDS ##########################################################
    ########################### BEADS, BACKBONEBEADS, AND MASS #############################################
    # Create Dictionary with Dataframes of Beads' Mass (easier for BBB lookup in calculation) and Bond List
    df_beadmass,bond_list,beads_list,backbonebeads_list = {}, [], [], []
    for l,lip in enumerate(LipidNames):
        if lip== 'POPC':
            Beads = [["C1A", "D2A", "C3A", "C4A", "GL1", "GL2", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"]]
            Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72]]           
            BackBoneBeads = [["PO4"]]
        elif lip == 'PUPC':
            Beads = [["D1A", "D2A", "D3A", "D4A", "D5A", "GL1", "GL2", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"]]        
            Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72]]
            BackBoneBeads = [["PO4"]]
        elif lip == 'POPE':
            Beads = [["C1A", "C2A", "C3A", "C4A", "GL1", "GL2", "C1B", "C2B", "D3B", "C4B", "C5B", "NH3", "PO4"]]
            Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72]]           
            BackBoneBeads = [["PO4"]]
        elif lip == 'PUPE':
            Beads = [["D1A", "D2A", "D3A", "D4A", "D5A", "GL1", "GL2", "C1B", "C2B", "C3B", "C4B", "NH3", "PO4"]]        
            Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72]]
            BackBoneBeads = [["PO4"]]
        elif lip == 'DPSM':
            Beads = [["AM2", "AM1", "T1A", "C2A", "C3A", "C1B", "C2B", "C3B", "C4B", "NC3", "PO4"]]        
            Mass = [[72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72]]
            BackBoneBeads = [["PO4"]]
        elif lip == 'CHOL':
            Beads = [["C2", "C1", "R5", "R4", "R3", "R2", "R1", "ROH"]]        
            Mass = [[72, 39.49, 0, 0, 159.65, 38.69, 0, 77.22]]
            BackBoneBeads = [["ROH"]]
        elif lip == 'ATOC':
            Beads = [["C4","C3","C2","C1","R3","R2","R1","ROH"]]
            Mass = [[72, 72, 72, 72, 72, 72, 72, 72]]
            BackBoneBeads = [["ROH"]]
        ######################################## BONDS ##########################################################
        phosphatidylcholine_bond_names = "NC3-GL1 "
        phosphatidylethanolamine_bond_names = "NH3-GL1 "
        # PCs
        if   lip == "POPC": bond_names = phosphatidylcholine_bond_names + "PO4-C1A PO4-GL2 GL1-D2A C1A-C3A D2A-C4A GL2-C2B C1B-C3B C2B-C4B GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
        elif lip == "PUPC": bond_names = phosphatidylcholine_bond_names + "PO4-D1A PO4-GL2 GL1-D2A D1A-D3A D2A-D4A D3A-D5A GL2-C2B C1B-C3B C2B-C4B GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
        # PEs
        elif lip == "PUPE": bond_names = phosphatidylethanolamine_bond_names + "PO4-D1A PO4-GL2 GL1-D2A D1A-D3A D2A-D4A D3A-D5A GL2-C2B C1B-C3B C2B-C4B GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
        elif lip == "POPE": bond_names = phosphatidylethanolamine_bond_names + "PO4-C1A PO4-GL2 GL1-C2A C1A-C3A C2A-C4A GL2-C2B C1B-D3B C2B-C4B D3B-C5B GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
        # DPSM
        elif lip == "DPSM": bond_names = "NC3-AM1 AM1-C2A AM2-C2B T1A-C3A C1B-C3B C2B-C4B AM2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
        # CHOL (newer V2_02 itp version)
        elif lip == "CHOL": bond_names = "ROH-C1 R3-C2 R2-C2\n"
        # ATOC
        # actual bonds of "ATOC": "ROH-R1 R1-R2 R2-R3 ROH-R3 R2-C1 ROH-R2 C1-C2 C2-C3 C3-C4"
        elif lip == "ATOC": bond_names = "ROH-C1 R1-C1 R3-C1 R2-C2 C1-C3 C2-C4 C1-C2 C2-C3 C3-C4"
        # Combine bond names
        bond_list = bond_list + [bond_names]
        # Dictionary with Dataframes of Beads' Mass
        d = {}
        d = {Beads[0][b]: [Mass[0][b]] for b in range(len(Beads[0]))}
        df_beadmass[lip] = pd.DataFrame(data=d)
        beads_list = beads_list+[Beads[0]]
        backbonebeads_list = backbonebeads_list + [BackBoneBeads[0]]
    return df_beadmass,bond_list,beads_list,backbonebeads_list
