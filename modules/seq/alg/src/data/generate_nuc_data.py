# code generation for subst_weight_matrix.cc
# loads NUC substitution matrix from file provided by NCBI at:
# ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/

def Generate(filename, matrix_name):

    with open(filename) as f:
        data = f.readlines()

    scores = dict()
    olcs = None    
    for line in data:
        if line.startswith('#'):
            continue
        if olcs is None:
            if " A " in line and " C " in line and " G " and " T " in line:
                # very high likelihood that this line contains the olcs
                olcs = line.strip().split()
            continue
        split_line = line.strip().split()
        if split_line[0] in olcs:
            olc = split_line[0]
            for score_idx, score in enumerate(split_line[1:]):
                scores[(olc, olcs[score_idx])] = score

    ost_olcs = ''.join(olcs)
    print("short %s[%i][%i]={"%(matrix_name, len(ost_olcs), len(ost_olcs)))
    for a in ost_olcs:
        score_string = "" 
        for b in ost_olcs:
            score = scores[(a,b)]
            if len(score_string) == 0:
                entry = [' ', ' ']
            else:
                entry = [',', ' ', ' ', ' ']
            for i in range(len(score)):
                entry[-1-i] = score[-1-i]
            score_string += ''.join(entry)
        print("  {%s},"%(score_string))
    print("};")
    print()

Generate("NUC.4.4_edited", "RAW_NUC44_DATA")

