import sys, os, glob
import csv


def write_csv(filename, darray):
    outfile = open(filename, 'w')
    wr = csv.writer(outfile)
    wr.writerows([['protein1', 'residue1', 'protein2', 'residue2', 'id', 'p value']])
    for entry in darray:
        if (entry[0] == '-') or (entry[0] == '*'):
            pass
        else:
            wr.writerows([entry])


    return True


def get_dist(r1, r2):
    import math
    x1, y1, z1 = r1[0], r1[1], r1[2]
    x2, y2, z2 = r2[0], r2[1], r2[2]
    return math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))

def parse_CAs(pdb_file):
    ca_array = {}
    with open(pdb_file) as fin:
        lines = fin.readlines()
    for line in lines:
        if line[12:16].strip() == 'CA':
            residue = int(line[22:26])
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            ca_array[residue] = [x, y, z]
    #print ca_array
    return ca_array

pdb_files = glob.glob("*.pdb")

print pdb_files
pvalue = 1.0
cutoffs = [ 5, 10, 15, 20, 25, 30]
delta = 0.15
for cutoff in cutoffs:
    count = 0
    print 'cutoff', cutoff
    xls = []
    for i in range(0,len(pdb_files)-1):
        pdb1 = pdb_files[i]
        for j in range(i+1, len(pdb_files)):
            pdb2 = pdb_files[j]
            print pdb1, pdb2
            pdb1_ca = parse_CAs(pdb1)
            pdb2_ca = parse_CAs(pdb2)
            pdb1_resi = pdb1_ca.keys()
            pdb2_resi = pdb2_ca.keys()

            for resi1 in pdb1_resi:
                coor1 = pdb1_ca[resi1]
                for resi2 in pdb2_resi:
                    coor2 = pdb2_ca[resi2]
                    dist = get_dist(coor1, coor2)
                    tpdb1 = pdb1.split("_")
                    tpdb2 = pdb2.split("_")
                    if (cutoff-delta) <= dist < (cutoff+delta):
                        count +=1
                        print tpdb1[0], resi1, tpdb2[0], resi2, count, pvalue
                        xls.append([tpdb1[0], resi1, tpdb2[0], resi2, count, pvalue])
                        #print xls

    write_csv('xl_'+str(cutoff)+".csv",xls)
