# code generation for subst_weight_matrix.cc
# loads BLOSUM substitution matrices from files provided by NCBI at:
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
            if " A " in line and " R " in line and " N " in line:
                # very high likelihood that this line contains the olcs
                olcs = line.strip().split()
            continue
        split_line = line.strip().split()
        if split_line[0] in olcs:
            olc = split_line[0]
            for score_idx, score in enumerate(split_line[1:]):
                scores[(olc, olcs[score_idx])] = score

    ost_olcs = "ABCDEFGHIKLMNPQRSTVWXYZ"
    print("short %s[23][23]={"%(matrix_name))
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

Generate("BLOSUM45", "RAW_BLOSUM45_DATA")
Generate("BLOSUM62", "RAW_BLOSUM62_DATA")
Generate("BLOSUM80", "RAW_BLOSUM80_DATA")
Generate("BLOSUM100", "RAW_BLOSUM100_DATA")

