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
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.mpi
import IMP.algebra
import IMP.pmi.analysis
import IMP.em
import IMP.core

# Import standard Python modules
import numpy as np
from math import log, pi, sqrt, exp
import math,csv
import collections
import re


#---------------------------
# Define Input Files
#---------------------------

datadirectory = "./"



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


