import sys, os, glob

pdb_file = sys.argv[1]
chain = sys.argv[2]


with open(pdb_file) as fin:
    lines = fin.readlines()


outfile = pdb_file.rstrip(".pdb")+"_"+chain+".pdb"

fout = open(outfile, 'w')

for line in lines:
    if line[0:4] == 'ATOM':
        pre_chain = line[0:21]
        chain = sys.argv[2]
        post_chain = line[22:]
        outline = pre_chain+chain+post_chain
        fout.write(outline)
    else:
        fout.write(line)


fout.close()
