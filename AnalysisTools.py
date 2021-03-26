
from os import listdir
import pandas as pd
import seaborn as sn
import numpy as np
import matplotlib.pyplot as plt

date_time = "20210321-224053"

def parseSmBits(date_time):
    kds = {'peptide86':0.7E-9,
           'peptide78':3.4E-9,
           'peptide79':8.5E-9,
           'peptide99':1.8E-7,
           'peptide128':2.8E-7,
           'native_test':0.9E-6,
           'peptide104':1.3E-6,
           'peptide101':2.5E-6,
           'peptide114':1.9E-4
           }

    l = listdir("./decoys/")
    fascs = []
    for f in l:
        if ".fasc" in f:
            if date_time in f:
                fascs.append(f)

    print(fascs)

    data = pd.read_json("./decoys/"+fascs[0], orient='records', lines=True)

    for f in fascs[1:]:
        d = pd.read_json("./decoys/"+f, orient='records', lines=True)
        data = data.append(d)

    print(data.shape)

    peptide_names = []
    for i, r in data.iterrows():
        s = r['filename'].split('-')
        s = s[0].split('/')
        name = s[2]
        peptide_names.append(name)

    data['peptide'] = peptide_names
    #print(data.head())

    peptides = set(data['peptide'])
    dfs = []

    for peptide in peptides:
        subdata = data[data['peptide']==peptide]
        subdata['kds'] = kds[peptide]

        stdev_total_score = np.std(subdata['total_score'])
        mean_total_score = np.mean(subdata['total_score'])
        print("{}: {} plus minus {}".format(peptide, mean_total_score, stdev_total_score))

        threshold_total_score = mean_total_score - 2 * stdev_total_score

        #dfs.append(subdata[subdata['total_score']<threshold_total_score])
        dfs.append(subdata.sort_values('total_score')[:10])

    significant_structures = pd.concat(dfs)
    print(significant_structures.shape)
    return significant_structures
