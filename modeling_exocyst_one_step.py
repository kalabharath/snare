import IMP
import RMF

import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf

import IMP.bayesianem
import IMP.bayesianem.restraint

import os, sys
import glob

import DLFCN as dl; 
sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)

import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof

import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking


#---------------------------
# Define Input Files
#---------------------------
datadirectory = "/netapp/sali/ilan/Exocyst/data/"
list_of_files = os.listdir("/netapp/sali/ilan/Exocyst/data/")
topology_file = datadirectory+"topology.txt"

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 10000
if '--test' in sys.argv: num_frames=20
num_mc_steps = 10

#--------------------------
# Create movers
#--------------------------
# rigid body movement params
rb_max_trans = 1.00
rb_max_rot = 0.01
# flexible bead movement
bead_max_trans = 2.00

#--------------------------------
# Build the Model Representation
#--------------------------------
# Initialize model
m = IMP.Model()

# Create list of components from topology file
beadsize = 10
colors = ['purple', 'blue', 'red', 'cornflower blue', 'light gray', 'turquoise', 'yellow', 'green']
components = ['Exo70', 'Exo84', 'Sec03', 'Sec05', 'Sec06', 'Sec08', 'Sec10', 'Sec15']
chains = "ABEFGHCD"
seqs = {}
mols = []
subs = []
# Create System and State
s = IMP.pmi.topology.System(m)
st = s.create_state()

# Setup sequences and chain ids quickly
offset=0

for n in range(len(components)):

    seqs = IMP.pmi.topology.Sequences("/netapp/sali/ilan/Exocyst/data/%s.txt" % components[n])
    mol = st.create_molecule(components[n], sequence=seqs['%s' % components[n]], chain_id=chains[n])

    for structure_file in list_of_files:
        if (structure_file.startswith(components[n]) and structure_file.endswith(".pdb")):
            
            print structure_file[0:-4]

            if structure_file == "Sec10_234_853.pdb":
                offset=5
            else:
                offset=0

            atomic = mol.add_structure("/netapp/sali/ilan/Exocyst/data/%s" % structure_file,  chain_id=chains[n], offset=offset)
            #print mol, list(atomic)[0], list(atomic)[-1]
            mol.add_representation(atomic, 
                                   resolutions=[1,10], 
                                   color=colors[n], 
                                   density_residues_per_component=10,
                                   density_prefix="/netapp/sali/ilan/Exocyst/data/gmm_structural_parts/"+structure_file[0:-4], 
                                   density_force_compute=False, 
                                   density_voxel_size=4.0)

    mols.append(mol)


'''
###Definition of the ideal helices
mols[0].add_representation(mols[0][10:35], resolutions=[1,10], ideal_helix=True, color=colors[0], density_residues_per_component=10,density_prefix="%s_gmm" % components[0], density_force_compute=False, density_voxel_size=4.0)
mols[1].add_representation(mols[1][230:270], resolutions=[1,10], ideal_helix=True, color=colors[1], density_residues_per_component=10,density_prefix="%s_gmm" % components[1], density_force_compute=False, density_voxel_size=4.0)
mols[1].add_representation(mols[1][315:330], resolutions=[1,10], ideal_helix=True, color=colors[1], density_residues_per_component=10,density_prefix="%s_gmm" % components[1], density_force_compute=False, density_voxel_size=4.0)
mols[3].add_representation(mols[3][945:960], resolutions=[1,10], ideal_helix=True, color=colors[3], density_residues_per_component=10,density_prefix="%s_gmm" % components[3], density_force_compute=False, density_voxel_size=4.0)
mols[4].add_representation(mols[4][355:369], resolutions=[1,10], ideal_helix=True, color=colors[4], density_residues_per_component=10,density_prefix="%s_gmm" % components[4], density_force_compute=False, density_voxel_size=4.0)
mols[5].add_representation(mols[5][76:95], resolutions=[1,10], ideal_helix=True, color=colors[5], density_residues_per_component=10,density_prefix="%s_gmm" % components[5], density_force_compute=False, density_voxel_size=4.0)
mols[5].add_representation(mols[5][100:130], resolutions=[1,10], ideal_helix=True, color=colors[5], density_residues_per_component=10,density_prefix="%s_gmm" % components[5], density_force_compute=False, density_voxel_size=4.0)
mols[6].add_representation(mols[6][74:115], resolutions=[1,10], ideal_helix=True, color=colors[6], density_residues_per_component=10,density_prefix="%s_gmm" % components[6], density_force_compute=False, density_voxel_size=4.0)
'''

###Definition of the beads regions
for mol in mols:
    color = ""
    if "Exo70" in mol.get_name():
        color = 'purple'
    elif "Exo84" in mol.get_name():
        color = 'blue'
    elif "Sec03" in mol.get_name():
        color = 'red'
    elif "Sec05" in mol.get_name():
        color = 'cornflower blue'
    elif "Sec06" in mol.get_name():
        color = 'light gray'
    elif "Sec08" in mol.get_name():
        color = 'turquoise'
    elif "Sec10" in mol.get_name():
        color = 'yellow'
    elif "Sec15" in mol.get_name():
        color = 'green'
    
    helices = []
    for hel in mol.get_ideal_helices():
        helices += hel
    mol.add_representation(mol[:]-mol.get_atomic_residues()-helices, resolutions=[10], color=color, setup_particles_as_densities=True)

'''
mols[0].add_representation(mols[0][:] - mols[0].get_atomic_residues()-mols[0][10:35], resolutions=[10], color=colors[0],setup_particles_as_densities=True)
mols[1].add_representation(mols[1][:] - mols[1].get_atomic_residues()-mols[1][200:275]-mols[1][315:330], resolutions=[10], color=colors[1],setup_particles_as_densities=True)
mols[2].add_representation(mols[2][:] - mols[2].get_atomic_residues()-mols[2][340:440]-mols[2][642:666]-mols[2][675:699]-mols[2][717:736], resolutions=[10], color=colors[2],setup_particles_as_densities=True)
mols[3].add_representation(mols[3][:] - mols[3].get_atomic_residues()-mols[3][460:490]-mols[3][528:551]-mols[3][945:960], resolutions=[10], color=colors[3],setup_particles_as_densities=True)
mols[4].add_representation(mols[4][:] - mols[4].get_atomic_residues()-mols[4][22:78]-mols[4][93:129]-mols[4][142:162]-mols[4][294:310]-mols[4][355:369], resolutions=[10], color=colors[4],setup_particles_as_densities=True)
mols[5].add_representation(mols[5][:] - mols[5].get_atomic_residues()-mols[5][76:95]-mols[5][100:130], resolutions=[10], color=colors[5],setup_particles_as_densities=True)
mols[6].add_representation(mols[6][:] - mols[6].get_atomic_residues()-mols[6][74:115], resolutions=[10], color=colors[6],setup_particles_as_densities=True)
mols[7].add_representation(mols[7][:] - mols[7].get_atomic_residues(), resolutions=[10], color=colors[7],setup_particles_as_densities=True)
'''
representation = s.build()
dof = IMP.pmi.dof.DegreesOfFreedom(m)

for mol in mols:
    print ""
    print "###ILAN###"
    print mol.get_name()


    selm = mol[:]-mol.get_atomic_residues()-helices
    dof.create_flexible_beads(selm,
                              max_trans=bead_max_trans)


    for structure_file in list_of_files:
        if structure_file.endswith(".pdb"):
            name=structure_file.strip("\n").split(".")[0].split("_")
            if name[0] in mol.get_name():
                print mol, name[0], name[1], name[2]
                sel0 = mol.residue_range(name[1], name[2])
                dof.create_rigid_body(sel0, 
                                      max_trans=rb_max_trans,
                                      max_rot=rb_max_rot,
                                      nonrigid_max_trans=bead_max_trans)
    dof.create_super_rigid_body(mol)


#--------------------------
# Define Degrees of Freedom
#--------------------------

# Add default mover parameters to simulation
outputobjects = [] # reporter objects (for stat files)
sampleobjects = [] # sampling objects

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files

# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan0", sf.evaluate(False)
EXO = []
crs = []

for mol in mols:
      cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol, scale=2.0)
      cr.add_to_model()
      outputobjects.append(cr)
      sampleobjects.append(cr)
      crs.append(cr)      
      EXO.append(mol)

print(EXO)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan0", sf.evaluate(False)

# Crosslinks - dataset
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging components together)
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("protein1")
kw.set_protein2_key("protein2")
kw.set_residue1_key("residue1")
kw.set_residue2_key("residue2")
kw.set_id_score_key(None)

xldb_exo = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb_exo.create_set_from_file(datadirectory+'Exocyst_MSXL.csv')

xls = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=representation,
                                                                            CrossLinkDataBase=xldb_exo,
                                                                            length=21,
                                                                            label="XLS",
                                                                            resolution=1.0,
                                                                            slope=0.02)


xls.rs.set_weight(25.0)
xls.add_to_model()
sampleobjects.append(xls)
outputobjects.append(xls)
dof.get_nuisances_from_restraint(xls)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan4", sf.evaluate(False)

IMP.pmi.tools.shuffle_configuration(representation)
dof.optimize_flexible_beads(200)

#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol

ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols, 
                                                              resolution=10)
ev1.add_to_model()
ev1.set_label('Exo')
outputobjects.append(ev1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan1", sf.evaluate(False)

densities = IMP.atom.Selection(representation,representation_type=IMP.atom.DENSITIES).get_selected_particles()
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,datadirectory+'map/model_g400.txt',slope=0.0000001,target_radii_scale=3.0,scale_target_to_mass=True)
gem.center_model_on_target_density(st)
gem.set_label("EM")
gem.add_to_model()
#sampleobjects.append(gem)
outputobjects.append(gem)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan1", sf.evaluate(False)


mc0=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=representation,
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
                                    replica_exchange_maximum_temperature=2.5,

                                    number_of_best_scoring_models=50,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=20000,
                                    global_output_directory='output/')


mc0.execute_macro()

rex0=mc0.get_replica_exchange_object()
