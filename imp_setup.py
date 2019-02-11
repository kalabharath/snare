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
all_cg_bead_size = 10
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

mols = []
subs = []


# Create system and state

s = IMP.pmi.topology.System(m)
st = s.create_state()

# Setup sequences and chain ids
offset = 0

print components

for n in range(0, len(components)):
    seqs = IMP.pmi.topology.Sequences(cwd+"/data/%s.fasta" %components[n])
    mol = st.create_molecule(components[n], sequence = seqs['%s' %components[n]], chain_id = chains[n])

    for pdb in list_of_files:
        if (pdb.startswith(components[n])) and pdb.endswith(".pdb"):
            atomic = mol.add_structure(cwd + '/data/%s' % pdb, chain_id=chains[n], offset=offset)
            mol.add_representation(atomic, resolutions=[1, 10], color=colors[n], bead_ca_centers=True)

    mols.append(mol)

for mol in mols:

    mol.add_representation(mol[:]-mol.get_atomic_residues(), resolutions = [all_cg_bead_size], color = colors[mols.index(mol)])


# calling System.build() creates all States and Molecules (and their representations)
# Once you call build(), anything without representation is destroyed.
representation = s.build()

# For verbose output of the representation
IMP.atom.show_with_representations(representation)


# --------------------------
# Define Degrees of Freedom
# --------------------------
# The DOF functions automatically select all resolutions
# Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
# Each structured unit is a single rigid body, with flexible beads corresponding to missing regions

dof = IMP.pmi.dof.DegreesOfFreedom(m)

for mol in mols:

    selm = mol[:] - mol.get_atomic_residues()
    dof.create_flexible_beads(selm, max_trans=bead_max_trans)

    for structure_file in list_of_files:
        if structure_file.endswith(".pdb"):
            name = structure_file.strip("\n").split(".")[0].split("_")
            if name[0] in mol.get_name():
                print mol, name[0], name[1], name[2]
                sel0 = mol.residue_range(name[1], name[2])
                dof.create_rigid_body(sel0,
                                      max_trans=rb_max_trans,
                                      max_rot=rb_max_rot,
                                      nonrigid_max_trans=bead_max_trans)

    # dof.create_super_rigid_body(mol)