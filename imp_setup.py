"""
This script samples SNARE complex proteins with fake cross-links data from known protein complexes
Author: Kala Bharath Pilla
Email: kalabharath@salilab.org
"""

import glob
import os
import sys

# Import IMP modules
import IMP.atom
import IMP.core
import IMP.pmi
import IMP.pmi.dof
import IMP.pmi.io.crosslink
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.tools
import IMP.pmi.topology

from vesicle_restraint import COMDistanceRestraint
from vesicle_restraint import VesicleMembraneRestraint

# ********************************************
#
# Define the SNARE complex system
#
# ********************************************

# Define Input Files

datadirectory = "./data/"
list_of_files = os.listdir("./data/")
all_cg_bead_size = 5
cwd = os.getcwd()
print cwd

# ---------------------------
# Define MC sampling Parameters
# ---------------------------

num_frames = 10000
if '--test' in sys.argv: num_frames = 20
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

# ********************************************
#
# Define the SNARE complex system
#
# ********************************************

colors = ['green', 'turquoise', 'magenta', 'orange']
components = ['Vamp2', 'Stx1a', 'Snap25', 'Cplx2']
chains = 'ABCD'

mols = []
subs = []

# Create system and state

s = IMP.pmi.topology.System(m, name="Snare")
st = s.create_state()

for n in range(0, len(components)):
    seqs = IMP.pmi.topology.Sequences(cwd + "/data/%s.fasta" % components[n])
    mol = st.create_molecule(components[n], sequence=seqs['%s' % components[n]], chain_id=chains[n])

    for pdb in list_of_files:
        if (pdb.startswith(components[n])) and pdb.endswith(".pdb"):
            atomic = mol.add_structure(cwd + '/data/%s' % pdb, chain_id=chains[n], offset=0)
            mol.add_representation(atomic, resolutions=[1, 5], color=colors[n], bead_ca_centers=True)
    mols.append(mol)

for mol in mols:
    mol.add_representation(mol[:] - mol.get_atomic_residues(), resolutions=[all_cg_bead_size],
                           color=colors[mols.index(mol)])

# calling System.build() creates all States and Molecules (and their representations)
# Once you call build(), anything without representation is destroyed.
h1_root = s.build()

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
                print "Defining DOFs", mol, name[0], name[1], name[2]
                sel0 = mol.residue_range(name[1], name[2])
                dof.create_rigid_body(sel0,
                                      max_trans=rb_max_trans,
                                      max_rot=rb_max_rot,
                                      nonrigid_max_trans=bead_max_trans)

# --------------------------
# Define Restraints
# --------------------------

# Here we are defining a number of restraints on our system.
# For all of them we call add_to_model() so they are incorporated into scoring
# We also add them to the outputobjects list, so they are reported in stat files

# Add default mover parameters to simulation
outputobjects = []  # reporter objects (for stat files)
sampleobjects = []  # sampling objects

# Excluded Volume Restraint
# To speed up this expensive restraint, we operate it at resolution 20

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))

snares = []

for mol in mols:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol, scale=2.0)
    cr.add_to_model()
    outputobjects.append(cr)
    sampleobjects.append(cr)
    snares.append(mol)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))

# --------------------------
# Membrane Restraints
# --------------------------


inside = [(266, 288, 'Stx1a'), (95, 116, 'Vamp2')]
above = [(1, 94, 'Vamp2'), (1, 265, 'Stx1a'), (1, 206, 'Snap25')]

mr = VesicleMembraneRestraint(h1_root, objects_inside=inside, objects_above=above, thickness=40, radius=100)
mr.add_to_model()
mr.create_membrane_density(file_out=cwd + "/vesicle.mrc")
outputobjects.append(mr)
dof.get_nuisances_from_restraint(mr)


# --------------------------
# Crosslinks - datasets
# --------------------------
# To use this restraint we have to first define the data format
# Here assuming that it's a CSV file with column names that may need to change
# Other options include the linker length and the slope (for nudging components together)

csv_files = glob.glob(cwd + "/data/*.csv")
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("protein1")
kw.set_protein2_key("protein2")
kw.set_residue1_key("residue1")
kw.set_residue2_key("residue2")
kw.set_id_score_key(None)

xls_objs = []

for csv_file in csv_files:
    xldb_exo = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
    xls_objs.append(xldb_exo)

for xldb in xls_objs:
    csv_file = csv_files[xls_objs.index(xldb)]
    xldb.create_set_from_file(csv_file)
    tlength = csv_file.split("/")
    tlength = tlength[-1].rstrip(".csv")
    tlength = int(tlength.lstrip("xl_"))

    xls = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=h1_root,
                                                                                CrossLinkDataBase=xldb,
                                                                                length=tlength,
                                                                                label="XLS_" + str(tlength),
                                                                                resolution=1.0,
                                                                                slope=0.02)
    xls.rs.set_weight(1.0)
    xls.add_to_model()
    sampleobjects.append(xls)
    outputobjects.append(xls)
    dof.get_nuisances_from_restraint(xls)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
sf.evaluate(False)


# ********************************************
#
# Define the Insulin vesicle system
#
# ********************************************

s2 = IMP.pmi.topology.System(m, name='Insulin Vesicle')
st2 = s2.create_state()
mol_v = st2.create_molecule('vesicle')
hier_bead = s2.build()

ch = hier_bead.get_children()
ves = IMP.atom.get_leaves(ch[0])[0]
xyzr = IMP.core.XYZR.setup_particle(ves.get_particle())
xyzr.set_coordinates_are_optimized(True)

xyzr.set_coordinates((0, 0, 400))
xyzr.set_radius(400)
IMP.atom.Mass.setup_particle(xyzr, 1)
IMP.atom.show_with_representations(hier_bead)

"""
# ------------------------------------
# Membrane restraint for the vesicle
# ------------------------------------
# probably not a good idea to implement: the vesicle Z coordinate is
# restrained and for the vesicle it would mean that it will be at the
# center of the sphere, unless you change the coordinate system
# *** this will not work so don't bother ***
vesicle_outputobjects = []
vesicle_mem_constraint = ['vesicle']
mr2 = VesicleMembraneRestraint(hier_bead, objects_inside=vesicle_mem_constraint, thickness=40, radius=100)
mr2.add_to_model()
print "wts inside mr2:", mr2
outputobjects.append(mr2)
dof.get_nuisances_from_restraint(mr2)
"""

# --------------------------
# Excluded Volume
# --------------------------
# This object defines all components to be sampled as well as the sampling protocol

ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols, other_objects=hier_bead,
                                                              resolution=10)
ev1.add_to_model()
ev1.set_label('Snare')
outputobjects.append(ev1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))

# -------------------------------------
# Combine hierarchies
# -------------------------------------

p = IMP.Particle(m)
hier_all = IMP.atom.Hierarchy.setup_particle(p)
hier_all.add_child(h1_root)
hier_all.add_child(hier_bead)
hier_all.set_name('System')

IMP.atom.show_with_representations(hier_all)

# -------------------------------------
# Distance restraint between
# -------------------------------------

print ves
print snares


vs_dist = COMDistanceRestraint(root_hier=hier_all, protein0='vesicle', protein1='Vamp2', distance=400, strength=1.0, label=None,
                               weight=1.0)
vs_dist.add_to_model()



# ---------------------------------
# Visualize initial configuration
# ---------------------------------

output = IMP.pmi.output.Output()
output.init_rmf("ini_all.rmf3", [hier_all])
output.write_rmf("ini_all.rmf3")
output.close_rmf("ini_all.rmf3")

print "Writing initial system state"
exit()
# ----------------------------------------------------
# Sampling the entire system by Monte-Carlo
# ----------------------------------------------------
# This object defines all components to be sampled as well as the sampling protocol

# IMP.pmi.tools.shuffle_configuration(h1_root)
# dof.optimize_flexible_beads(200)

mc0 = IMP.pmi.macros.ReplicaExchange0(m,
                                      root_hier=hier_all,
                                      monte_carlo_sample_objects=dof.get_movers(),
                                      output_objects=outputobjects,
                                      crosslink_restraints=sampleobjects,
                                      monte_carlo_temperature=1.0,

                                      simulated_annealing=True,
                                      simulated_annealing_minimum_temperature=1.0,
                                      simulated_annealing_maximum_temperature=1.5,
                                      simulated_annealing_minimum_temperature_nframes=200,
                                      simulated_annealing_maximum_temperature_nframes=20,

                                      replica_exchange_minimum_temperature=1.0,
                                      replica_exchange_maximum_temperature=5.0,

                                      monte_carlo_steps=num_mc_steps,
                                      number_of_frames=20000,
                                      global_output_directory='output/')

mc0.execute_macro()

rex0 = mc0.get_replica_exchange_object()
