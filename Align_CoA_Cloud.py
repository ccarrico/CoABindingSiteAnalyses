#! /usr/bin/python

# Crude parallelizer: sys.argv = [this.py] [K_N] where k is this process's ID (0-7) and N is number of processes (8)

import sys
import time
import os
sys.path.append('/home/ccarrico/pcmd/')
from AtomTools import *

print "Run with arguments:", sys.argv
[Kprocess, Nprocess] = [int(x) for x in sys.argv[1].split('_')]
processCounter = 0

pdb_filename    = 'Q9WUM5_1EUD/Q9WUM5_1EUD.pdb'

ligand_filename = 'Q9WUM5_1EUD/2SCU_D_301_COA.pdb'

CoA_cloud_fname = 'CoA_cloud.pdb'

pdb_set = CoordSet([])
pdb_set.pdb_load(pdb_filename)

ligand_set = CoordSet([])
ligand_set.pdb_load(ligand_filename)

CoA_set = CoordSet([])
CoA_set.pdb_load(CoA_cloud_fname)
CoA_residue_numbers = [x for x in set([x.residue_number for x in CoA_set.Coords])] #list copy so it can be repeatably indexed

CoA_Alignment_Trio = [CoA_set.get_res_number(1).get_atom_type('P1A').Coords[0],
                      CoA_set.get_res_number(1).get_atom_type('O3A').Coords[0],
                      CoA_set.get_res_number(1).get_atom_type('P2A').Coords[0]]

filename_list_lines = file('Good_COA_RMS_placements.txt').read().split("\n")[:-1]
LigandDict = {}

for list_line in filename_list_lines:

    [fname, ligand, rms] = list_line.split()
    if not fname in LigandDict:
        LigandDict.update({fname:[]})
    LigandDict[fname].append(ligand)

GlobalTupleStericDict = {}
GlobalTupleClashSpacing = max(VDW_Dictionary.values())

MAX_CLASH = 2.0

StartTime = time.time()

for filename in LigandDict:

    pdb_set.pdb_load(filename + '/' + filename + '.pdb')
    folder = filename

    CoA_ligand_list = [filename + '/' + x for x in LigandDict[filename]]

    print filename, len(CoA_ligand_list)

    for ligand_fname in CoA_ligand_list:

        if os.path.isfile(ligand_fname.replace('.pdb','_CoA_cloud.pdb')):
            continue

        processCounter +=1 
        if processCounter % Nprocess != Kprocess:
            continue
        else:
            print 'Running on', ligand_fname, 'for', filename
        
        ligand_set.pdb_load(ligand_fname)

        try:
            Ligand_Alignment_Trio = [ligand_set.get_atom_type('P1A').Coords[0],
                                     ligand_set.get_atom_type('O3A').Coords[0],
                                     ligand_set.get_atom_type('P2A').Coords[0]]
        except IndexError:
            continue

        CoA_set.align(CoA_Alignment_Trio, Ligand_Alignment_Trio)

        Valid_CoA_resnums = [x for x in CoA_residue_numbers if CoA_set.get_res_number(x).clash_check(pdb_set) < MAX_CLASH]
        CoA_set.get_res_numbers(Valid_CoA_resnums).dump(ligand_fname.replace('.pdb','_CoA_cloud.pdb'))

        print filename, ligand_fname, len(Valid_CoA_resnums), time.time() - StartTime
        
print "done"

