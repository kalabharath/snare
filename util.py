def read_struc_file(fn):

    """Read distances predicted from NOEs as NOEs from file."""
    coupRes = []

    with open(fn, "rt") as f:
        for i, row in enumerate(csv.reader(f, delimiter="\t")):
            # interactList.append([row[0].split(":")[1], row[1], row[2], row[3], row[9]])
            print(row)
            coupRes.append([int(row[0]), int(row[2])])
    return coupRes
