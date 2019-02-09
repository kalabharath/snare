"""
This script samples SNARE complex proteins with fake cross-links data from known complexes
Author: Kala Bharath Pilla
Email: kalabharath@salilab.org

"""

# Import IMP modules
import IMP
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.crosslinking
import IMP.algebra
import IMP.pmi.analysis

import IMP.core

# Import standard Python modules
import numpy as np
from math import log, pi, sqrt, exp
import math,csv
import collections
import re
import sys, os, glob


# ---------------------------
# Define Input Files
# ---------------------------

datadirectory = "./data/"
list_of_files = os.listdir("./data/")
topology_file = datadirectory+"topology.txt"
cwd = os.getcwd()
print cwd
# ---------------------------
# Define MC sampling Parameters
# ---------------------------

num_frames = 10000
if '--test' in sys.argv: num_frames=20
num_mc_steps = 10


# --------------------------
# Create movers
# --------------------------
# rigid body movement params
rb_max_trans = 1.00
rb_max_rot = 0.01
# flexible bead movement
bead_max_trans = 2.00

# --------------------------------
# Build the Model Representation
# --------------------------------
# Initialize model
m = IMP.Model()


# create list of components from topology file

beadsize = 10
colors = ['green', 'turquoise', 'pink', 'yellow', 'red']
components = ['Vamp2', 'Stx1a', 'Snap25', 'Cplx2', 'Syt7']
chains = 'ABCDE'

seqs = {}
mols = {}
subs = []


# Create system and state

s = IMP.pmi.topology.System(m)
st = s.create_state()



# Setup sequences and chain ids

offset = 0


for n in range(0, len(components)):
    seqs = IMP.pmi.topology.Sequences(cwd+"/data/%s.fasta" %components[n])
    mol = st.create_molecule(components[n], sequence = seqs['%s' %components[n]], chain_id = chains[n])

    for pdb in list_of_files:
        if (pdb.startswith(components[n])) and pdb.endswith(".pdb"):
            print pdb, 'wtf'
            atomic = mol.add_structure(cwd+'/data/%s' %pdb, chain_id = chains[n], offset= offset)


"""


# Input sequences
spots_seq_file = '/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/spots.fasta'
spots_seqs = IMP.pmi.topology.Sequences(spots_seq_file)
spots_components={"LCB1":["A","F"],"LCB2":["B","G"],"ORM1":["C","H"],"ORM2":["D","I"],"TSC3":["E","J"]}
spots_colors ={"LCB1":["blue","cyan"],"LCB2":["red","salmon"],"ORM1":["green","gold"],"ORM2":["pink","orange"],"TSC3":["brown","black"]}


# input files
allcrosslink_file = "/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/crossLinks/allCrossLinks.csv"

# Parameters to tweak for production run
all_cg_bead_size = 10

# MC parameters
RB_MAX_TRANS = 4.0
RB_MAX_ROT = 1.0
FLEX_MAX_TRANS = 4.0
SRB_MAX_TRANS = 1.0
SRB_MAX_ROT = 0.1


"""