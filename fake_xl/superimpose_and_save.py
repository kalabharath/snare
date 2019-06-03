import sys, os, glob
import pymol
pymol.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
pymol.finish_launching()

pdb_files = glob.glob("*_*_*.pdb")
print (pdb_files)

for pdb in pdb_files:
    native = pymol.cmd.load("complex.cif", 'obj1')
    pymol.cmd.load(pdb, 'obj2')
    pymol.cmd.align('obj2', 'obj1')
    pymol.cmd.save(pdb, 'obj2', format='pdb')
    pymol.cmd.reinitialize()
