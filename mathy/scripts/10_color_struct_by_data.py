#!/usr/bin/env python
# this script should be in the scripts/ subdirectory of the project

'''This script colors the Ran structure from PDB ID 1K5D according to the
number of beneficial or deleterious mutations at each position. A mutation is
considered beneficial if the score in the 15 generation EMPIRIC dataset
is > 1, and deleterious if it is < -1.

The number of mutations in each category are counted up for each position, and
the color scale is from 0-10, so a position with 11 deleterious mutations is
just as red as a position with 15 deleterious mutations.

The coloring is done separately for beneficial and deleterious mutations, so a
position that has 5 deleterious mutations and 4 beneficial ones will be colored
for the 5 in the red structure and the 4 in the blue structure.

adapted from:
http://betainverse.github.io/blog/2014/10/13/pymol-color-by-data/
'''

# import packages
import os, sys 
import pandas as pd, numpy as np
import pymol
import pyrosetta
from pyrosetta import *
from pyrosetta_utils import aminoacids


# First, we prepare a file with scores to be read in by the other PyMOL functions
# We need this file to be tab-delimited with pdb residues and scores:
# 4   0.525
# 5   0.335
# 6   0.356

# read in the pdb to be colored using pyrosetta, so we can know which residues
# are present in the structure
pyrosetta.init()
pdb = '../data/pdbs_ran/1k5d.pdb'
pose = pyrosetta.pose_from_pdb(pdb)
resnums_pose = [res.seqpos() for res in pose.residues if res.name1() in aminoacids.ONELETTER]
resnums_pdb = [int(pose.pdb_info().pose2pdb(n).split(' ')[0]) for n in resnums_pose]

# name datafiles
bene_data = 'b_temp.txt'
dele_data = 'd_temp.txt'

# make a dataframe that counts the number of beneficial and deleterious muts
df = (pd.read_csv('../data/empiric_data_cleaned.txt', sep='\t')
        .query('dataset == "15gen"')
        .assign(is_beneficial = lambda df: df.fitness_score.apply(
            lambda x: True if x > 1 else False),
                is_deleterious = lambda df: df.fitness_score.apply(
            lambda x: True if x < -1 else False))
        .set_index('position')
        .loc[resnums_pdb][['is_beneficial','is_deleterious']]
        .groupby('position')
        .agg(np.sum)
        .apply(pd.to_numeric, downcast = 'integer')
     )

# make datafiles
df[['is_beneficial']].to_csv(bene_data, sep='\t', header=False)
df[['is_deleterious']].to_csv(dele_data, sep='\t', header=False)

# load in pdb, and set view
def load_into_pymol(name):
    pymol.cmd.load(pdb, name)
    pymol.cmd.set_view('''
         0.141285628,    0.062554061,   -0.987988949,\
        -0.848635197,    0.521547079,   -0.088338152,\
         0.509757161,    0.850924432,    0.126771584,\
        -0.000294507,   -0.000136614, -224.959335327,\
        62.194507599,    1.145602703,   84.127868652,\
       177.360321045,  272.560058594,  -20.000000000''')
    pymol.cmd.set('valence', 'off')
    pymol.cmd.set('depth_cue', 'off')
    pymol.cmd.set('specular', 'off')
    pymol.cmd.set('ray_shadows', 'off')
    pymol.cmd.set('surface_quality', '3')
    pymol.cmd.hide('everything', name+' and resn HOH')
    pymol.cmd.color('atomic', name+' and not polymer')
    pymol.cmd.color('salmon', name+' and not polymer and elem C')
    return

load_into_pymol('n_beneficial_mutations')
load_into_pymol('n_deleterious_mutations')

# load in pymol scripts for coloring based on b factor
pymol.cmd.run('pymol_scripts/data2bfactor.py')
pymol.cmd.run('pymol_scripts/spectrumany.py')

# prepare output image names
bene_out = '../plots/struct_colored_beneficial'
dele_out = '../plots/struct_colored_deleterious'

# color beneficial 
pymol.cmd.hide('everything', 'n_deleterious_mutations')
pymol.cmd.select('residues', 'n_beneficial_mutations and polymer')
pymol.cmd.alter('residues', 'b=0')
pymol.cmd.do('data2b_res n_beneficial_mutations, ' + bene_data)
pymol.cmd.spectrum('b', 'white blue', 'residues', 0, 10)
pymol.cmd.png(bene_out+'.png')
pymol.cmd.rotate('y', '180')
pymol.cmd.png(bene_out+'_180.png')

# color deleterious
pymol.cmd.rotate('y', '180')
pymol.cmd.hide('everything', 'n_beneficial_mutations')
pymol.cmd.show('cartoon', 'n_deleterious_mutations')
pymol.cmd.show('sticks', 'n_deleterious_mutations and not polymer')
pymol.cmd.show('nb_spheres', 'n_deleterious_mutations and resn Mg')
pymol.cmd.select('residues', 'n_deleterious_mutations and polymer')
pymol.cmd.alter('residues', 'b=0')
pymol.cmd.do('data2b_res n_deleterious_mutations, ' + dele_data)
pymol.cmd.spectrum('b', 'white red', 'residues', 0, 10)
pymol.cmd.png(dele_out+'.png')
pymol.cmd.rotate('y', '180')
pymol.cmd.png(dele_out+'_180.png')

# save session file
pymol.cmd.show('cartoon', 'n_beneficial_mutations')
pymol.cmd.show('sticks', 'n_beneficial_mutations and not polymer')
pymol.cmd.show('nb_spheres', 'n_beneficial_mutations and resn Mg')
pymol.cmd.set('grid_mode', value=1)
pymol.cmd.deselect()
pymol.cmd.save('../plots/mutation_counts_pymol.pse')


# delete the datafiles since they're not well annotated anyway
os.remove(bene_data)
os.remove(dele_data)




