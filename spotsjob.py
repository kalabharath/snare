"""This script samples SPOTS complex subunits with cross-links data and Cryo-EM 3D map from Adam Frost's group"""

import IMP
import RMF
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
import IMP.pmi.restraints.em
import IMP.pmi.restraints.crosslinking
import IMP.bayesianem
import IMP.bayesianem.restraint
import IMP.mpi
import os,sys
import IMP.algebra
import tempfile
import numpy as np
import IMP.pmi.analysis
import IMP.em
import IMP.core
import sys
import os
from math import log, pi, sqrt, exp
#from LognormalAmbiguousRestraint import autoLognormalAmbiguousRestraint
import math,csv
import collections
import re


def read_struc_file(fn):
    """Read distances predicted from NOEs as NOEs from file."""
    coupRes = []
    with open(fn, "rt") as f:
        for i, row in enumerate(csv.reader(f, delimiter="\t")):
            #interactList.append([row[0].split(":")[1], row[1], row[2], row[3], row[9]])
            print(row)
            coupRes.append([int(row[0]), int(row[2])])
    return coupRes


# Function for creating representation of subunits without any structure

def add_spots_rep(mol, chain, unstructured_bead_size, clr, prot, dens):

    print(mol, chain, unstructured_bead_size, clr, prot, dens)

    if prot=='LCB1' or prot=='LCB2':

        atomic = mol.add_structure('/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/structures/LCB1_LCB21.pdb', chain_id=chain, offset=0)

        mol.add_representation(atomic, resolutions=[1,10], color = clr, density_prefix='/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/em_data/'+ prot +'.'+ chain, density_residues_per_component=dens, density_force_compute=False)

        mol.add_representation(mol[:] - atomic, resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

    # if prot=='LCB2':

    #   atomic = mol.add_structure('/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/structures/LCB2.pdb', chain_id=chain, offset=0)

    #   mol.add_representation(atomic, resolutions=[1,10], color = clr, density_prefix='/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/em_data/'+ prot +'.'+ chain, density_residues_per_component=dens, density_force_compute=False)

    #   mol.add_representation(mol[:] - atomic, resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

    if prot=='ORM1':

        atomic = mol.add_structure('/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/structures/'+ prot +".pdb", chain_id=chain, offset=0)

        mol.add_representation(atomic, resolutions=[1,10], color = clr, density_prefix='/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/em_data/'+ prot +'.'+ chain, density_residues_per_component=dens, density_force_compute=False)

        mol.add_representation(mol[:] - atomic, resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

    if prot=='ORM2':

        mol.add_representation(mol[:81],resolutions=[unstructured_bead_size],color=clr,setup_particles_as_densities=True,density_force_compute=False)

        mol.add_representation(mol[81:95], resolutions=[1,3], color=clr, density_force_compute=False, ideal_helix=True)

        mol.add_representation(mol[95:109], resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

        mol.add_representation(mol[109:122], resolutions=[1,3], color=clr, density_force_compute=False, ideal_helix=True)

        mol.add_representation(mol[122:158], resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

        mol.add_representation(mol[158:167], resolutions=[1,3], color=clr, density_force_compute=False, ideal_helix=True)

        mol.add_representation(mol[167:178], resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

        mol.add_representation(mol[178:190], resolutions=[1,3], color=clr, density_force_compute=False, ideal_helix=True)

        mol.add_representation(mol[190:], resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

    if prot=='TSC3':

        atomic = mol.add_structure('/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/structures/'+ prot +".pdb", chain_id=chain, offset=0)

        mol.add_representation(atomic, resolutions=[1,10], color = clr, density_prefix='/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/em_data/'+ prot +'.'+ chain, density_residues_per_component=dens, density_force_compute=False)

        mol.add_representation(mol[:59] - atomic, resolutions=[unstructured_bead_size], color=clr, setup_particles_as_densities=True, density_force_compute=False)

        mol.add_representation(mol[59:80], resolutions=[1,3], color=clr, density_force_compute=False, ideal_helix=True)

    return mol

# Function for creating rotational symmetry
def create_rotational_symmetry(mdl, mol1, mol2, transform,resolution='all'):
    '''
    The copies must not contain rigid bodies.
    The symmetry restraints are applied at each leaf.
    '''
    sm = IMP.core.TransformationSymmetry(transform)

    href    = IMP.pmi.tools.input_adaptor(mol1,resolution,flatten=True)
    hclones = IMP.pmi.tools.input_adaptor(mol2,resolution,flatten=True)

    # Remove densities
    hhref    = [h for h in href if (len(IMP.atom.Fragment(h).get_residue_indexes())>0 or 'Residue' in h.get_name() or 'bead' in h.get_name())]
    hhclones = [h for h in hclones if (len(IMP.atom.Fragment(h).get_residue_indexes())>0 or 'Residue' in h.get_name() or 'bead' in h.get_name())]

    ref_rbs,ref_beads = IMP.pmi.tools.get_rbs_and_beads(hhref)
    clones_rbs,clones_beads = IMP.pmi.tools.get_rbs_and_beads(hhclones)

    # Check that the clones do not contain any rigid bodies
    if len(clones_rbs)>0:
        raise Exception("ERROR: Clones should not contain rigid bodies")
    if len(hhref)!=len(hhclones):
        raise Exception("ERROR: Your references don't match your clones")

    lc = IMP.container.ListSingletonContainer(mdl)

    for n, p in enumerate(hhref):
        pc = hhclones[n]
        IMP.core.Reference.setup_particle(pc.get_particle(), p.get_particle())
        lc.add(pc.get_particle().get_index())

    c = IMP.container.SingletonsConstraint(sm, None, lc)
    mdl.add_score_state(c)

    dof.disable_movers(hclones)

###################### SYSTEM SETUP #####################
# Parameters to tweak for production run

#input files
allcrosslink_file = "/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/crossLinks/allCrossLinks.csv"

all_cg_bead_size = 10

#MC parameters
RB_MAX_TRANS = 4.0
RB_MAX_ROT = 1.0
FLEX_MAX_TRANS = 4.0
SRB_MAX_TRANS = 1.0
SRB_MAX_ROT = 0.1

# Input sequences
spots_seq_file = '/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/spots.fasta'
spots_seqs = IMP.pmi.topology.Sequences(spots_seq_file)
spots_components={"LCB1":["A","F"],"LCB2":["B","G"],"ORM1":["C","H"],"ORM2":["D","I"],"TSC3":["E","J"]}
spots_colors ={"LCB1":["blue","cyan"],"LCB2":["red","salmon"],"ORM1":["green","gold"],"ORM2":["pink","orange"],"TSC3":["brown","black"]}

# Setup System and add a State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

# Add Molecules for each component as well as representations
mols = []
spots_mols = []

# Loop over to create two copies for subunits with partial structural coverage
j=0
for prot in ['LCB1', 'LCB2', 'ORM1', 'ORM2', 'TSC3']:
  for i,chain in enumerate(spots_components[prot]):
    if i==0:
      mol1 = st.create_molecule(prot, sequence=spots_seqs[prot + ":" + chain], chain_id=chain)
    else:
      mol1 = spots_mols[j].create_copy(chain_id=chain)

    color = spots_colors[prot][i]
    mol = add_spots_rep(mol1, chain, all_cg_bead_size, color, prot, 10)
    mols.append(mol)
    spots_mols.append(mol)
  j = j + 2

# calling System.build
# calling System.build() creates all States and Molecules (and their representations)
# Once you call build(), anything without representation is destroyed.
# You can still use handles like molecule[a:b], molecule.get_atomic_residues() or molecule.get_non_atomic_residues()
# However these functions will only return BUILT representations
root_hier = s.build()

# For verbose output of the representation
IMP.atom.show_with_representations(root_hier)

##############################################################
# Setup degrees of freedom
# The DOF functions automatically select all resolutions
# Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

# DOF for SPOTS complex proteins
# Each structured unit is a single rigid body, with flexible beads corresponding to missing regions

lcb_dimer_u = [spots_mols[0].get_non_atomic_residues(),  spots_mols[2].get_non_atomic_residues()]
lcb_dimer = [spots_mols[0],  spots_mols[2]]
dof.create_rigid_body(lcb_dimer, nonrigid_parts=lcb_dimer_u, max_trans=RB_MAX_TRANS, max_rot=RB_MAX_ROT, nonrigid_max_trans=FLEX_MAX_TRANS, name="spots_"+ str(i) +"_rb")

for i, mol in enumerate(spots_mols):
# Create a rigid body for each domain with structural information, the flexible beads inside are part of the rigid body, the remaining part of the protein is a flexible object

# LCB1 has only partial structure from comparative modeling
    # if i==0:

    #   dof.create_rigid_body(mol, nonrigid_parts=mol.get_non_atomic_residues(), max_trans=RB_MAX_TRANS, max_rot=RB_MAX_ROT, nonrigid_max_trans=FLEX_MAX_TRANS, name="spots_"+ str(i) +"_rb")

    # if i==2:

    #   dof.create_rigid_body(mol, nonrigid_parts=mol.get_non_atomic_residues(), max_trans=RB_MAX_TRANS, max_rot=RB_MAX_ROT, nonrigid_max_trans=FLEX_MAX_TRANS, name="spots_"+ str(i) +"_rb")

    if i==4:

        dof.create_rigid_body(mol, nonrigid_parts=mol.get_non_atomic_residues(), max_trans=RB_MAX_TRANS, max_rot=RB_MAX_ROT, nonrigid_max_trans=FLEX_MAX_TRANS, name="spots_"+ str(i) +"_rb")

    # if i==6:

    #   dof.create_rigid_body(mol, nonrigid_parts=mol.get_non_atomic_residues(), max_trans=RB_MAX_TRANS, max_rot=RB_MAX_ROT, nonrigid_max_trans=FLEX_MAX_TRANS, name="spots_"+ str(i) +"_rb")

    if i==6:

        dof.create_rigid_body(mol[81:95], max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,nonrigid_max_trans=FLEX_MAX_TRANS,name="spots_"+ str(i) +"_rb")

        dof.create_rigid_body(mol[109:122], max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,nonrigid_max_trans=FLEX_MAX_TRANS,name="spots_"+ str(i) +"_rb")

        dof.create_rigid_body(mol[158:167], max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,nonrigid_max_trans=FLEX_MAX_TRANS,name="spots_"+ str(i) +"_rb")

        dof.create_rigid_body(mol[178:190], max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,nonrigid_max_trans=FLEX_MAX_TRANS,name="spots_"+ str(i) +"_rb")

        dof.create_rigid_body(mol.get_non_atomic_residues(), nonrigid_parts=mol.get_non_atomic_residues(), max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,nonrigid_max_trans=FLEX_MAX_TRANS,name="spots_"+ str(i) +"_rb")

    if i==8:

        dof.create_rigid_body(mol[1:59], nonrigid_parts=mol.get_non_atomic_residues(), max_trans=RB_MAX_TRANS, max_rot=RB_MAX_ROT, nonrigid_max_trans=FLEX_MAX_TRANS, name="spots_"+ str(i) +"_rb")

        dof.create_rigid_body(mol[59:80], max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,nonrigid_max_trans=FLEX_MAX_TRANS,name="spotshelix_"+ str(i) +"_rb")

####################### RESTRAINTS #####################
output_objects = [] # keep a list of functions that need to be reported
display_restraints = [] # display as springs in RMF

# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
for mol in mols:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.add_to_model()
    output_objects.append(cr)
    crs.append(cr)

# Excluded volume - automatically more efficient due to rigid bodies. The excluded volume here is evaluated
mols_sim = [spots_mols[0], spots_mols[2], spots_mols[4], spots_mols[6], spots_mols[8]]
mols_sym = [spots_mols[1], spots_mols[3], spots_mols[5],spots_mols[7], spots_mols[9]]

evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols_sim)
evr1.add_to_model()
evr1.set_label('Intra')
output_objects.append(evr1)

evr2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols_sim, other_objects=mols_sym)
evr2.add_to_model()
evr2.set_label('Inter')
output_objects.append(evr2)

###############################
# Symmetry constraints

# C2 symmetry
center = IMP.algebra.Vector3D([0,0,0])
ref = IMP.algebra.get_rotation_about_axis([0.0,0.0,1.0],pi)

# Apply the C2 symmetry transformation and moves along this symmetry axis
transform_ref = IMP.algebra.get_rotation_about_point(center,ref)
create_rotational_symmetry(mdl,spots_mols[0], spots_mols[1],transform_ref)
create_rotational_symmetry(mdl,spots_mols[2], spots_mols[3],transform_ref)
create_rotational_symmetry(mdl,spots_mols[4], spots_mols[5],transform_ref)
create_rotational_symmetry(mdl,spots_mols[6], spots_mols[7],transform_ref)
create_rotational_symmetry(mdl,spots_mols[8], spots_mols[9],transform_ref)


###############################
# Membrane Restraint

mr = IMP.pmi.restraints.basic.MembraneRestraint(root_hier,
                                                         objects_inside=[(55, 75, 'LCB1'), (59, 79, 'LCB2'), (84, 104, 'ORM1'), (110, 130, 'ORM1'), (163, 183, 'ORM1'), (185, 205, 'ORM1'), (78, 98, 'ORM2'), (104, 124, 'ORM2'), (157, 177, 'ORM2'), (179, 199, 'ORM2'), (59, 80, 'TSC3')],
                                                         objects_above=[(1, 54, 'LCB1'), (76, 552, 'LCB1'), (1, 58, 'LCB2'), (80, 524, 'LCB2'), (59, 79, 'ORM1'), (1, 83, 'ORM1'), (105, 109, 'ORM1'), (131, 162, 'ORM1'), (206, 222, 'ORM1'), (1, 77, 'ORM2'), (99, 103, 'ORM2'), (125, 156, 'ORM2'), (200, 216, 'ORM2'), (1, 58, 'TSC3')],thickness=40)


mr.add_to_model()
output_objects.append(mr)
dof.get_nuisances_from_restraint(mr)

###############################
#Crosslink restraint
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")

# Intra-molecular crosslinks
allcrosslinkdb = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
allcrosslinkdb.create_set_from_file(allcrosslink_file)
allcrosslinkxlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hier,CrossLinkDataBase=allcrosslinkdb,length=21.0,label="XLDSS",filelabel='dss',resolution=1,slope=0.01)
allcrosslinkxlr.add_to_model()
allcrosslinkxlr.set_weight(75)
output_objects.append(allcrosslinkxlr)
display_restraints.append(allcrosslinkxlr)

# Sampling of All crosslinks psi values
allcrosslinkxlr.set_psi_is_sampled(True)
psi = allcrosslinkxlr.psi_dictionary["PSI"][0]

dof.get_nuisances_from_restraint(allcrosslinkxlr) # needed to sample the nuisance particles (noise params)

################################
# Coevolutionary Restraint
#tuple_selection_pairs = read_struc_file("/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/coevol/TSC3.txt")

#rnoes = autoLognormalAmbiguousRestraint(tuple_selection_pairs, root_hier)
#rnoes.add_to_model()
#output_objects.append(rnoes)
#dof.get_nuisances_from_restraint(rnoes)

#################################
#cryo-EM restraints
densities = IMP.atom.Selection(root_hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,'/wynton/scratch/rakesh/SPOTSModeling/spotsJobLarge/data/Centered_Map_Cropped.txt',slope=0.000001,target_radii_scale=3.0,scale_target_to_mass=True)
gem.set_label("EM")
gem.add_to_model()
output_objects.append(gem)

####################### SAMPLING #####################
for state in st.system.get_states():
    IMP.pmi.tools.shuffle_configuration(state.get_hierarchy(),
                                        max_translation=300)
dof.optimize_flexible_beads(100)

mc = IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    crosslink_restraints=display_restraints,                     # will display like XLs
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0 ,
                                    num_sample_rounds = 1,
                                    number_of_best_scoring_models=0,
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='jobResults/',
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=100000
                                   )

mc.execute_macro()
