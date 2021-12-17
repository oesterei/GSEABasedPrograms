# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:54:06 2021
@author: Laura K. Harris, Ph.D.
"""

import numpy as np
import pandas as pd
import scipy
import gseapy as gp
from gseapy.plot import gseaplot

def zscore(df, datasetnum):
    tempnum = 0
    df2 = pd.DataFrame()
    totalrows = df[df.columns[0]].count()

    while tempnum < (totalrows):
        zscorearray = []
        temparray = df.iloc[tempnum].to_numpy()
        geneID = temparray[0]
        temparray = np.delete(temparray, 0, 0)

        mean = np.mean(temparray)
        std = np.std(temparray)

        for value in temparray:
            tempz = (value - mean) / std
            zscorearray = np.append(zscorearray, tempz)
        tempseries = pd.Series(zscorearray, name = geneID)
        df2 = df2.append(tempseries)
        
        if tempnum % 1000 == 0:
            print('Z-scoring gene #' + str(tempnum))
        tempnum = tempnum + 1
    
    df2copy = df2
    header = list(df.head(0))
    header.pop(0)
    df2.columns = header
    
    filename = 'Zscoredata' + str(datasetnum) + '.txt'
    df2.to_csv(filename, sep='\t', index=True) #change file name as needed
    print('Z-score done for dataset #' + str(datasetnum))
    print()
    return df2copy

def Tscore(df, controlcols, experimentalcols, datasetnum):
    dfstat = pd.DataFrame()
    dfpval = pd.DataFrame()
    tempnum = 0
    for r in df.iterrows():
        geneID = r[0]
        row = r[1]
        controlvals = []
        experimentalvals = []
        for iter in controlcols:
            controlvals.append(row[iter])
        for iter in experimentalcols:
            experimentalvals.append(row[iter])
        stat, pval = scipy.stats.ttest_ind(experimentalvals, controlvals, equal_var=False)
        tempstatseries = pd.Series(stat, name = geneID)
        temppvalseries = pd.Series(pval, name = geneID)
        dfstat = dfstat.append(tempstatseries)
        dfpval = dfpval.append(temppvalseries)
        
        if tempnum % 1000 == 0:
            print('T-scoring gene #' + str(tempnum))
        tempnum = tempnum + 1
        
    dfstat.columns = ['Tscore']
    df2 = pd.concat([dfstat, dfpval], axis = 1)
    df2.columns = ['Tscore', 'pval']
    filename = 'Tscoredata' + str(datasetnum) + '.txt'
    df2.to_csv(filename, sep='\t', index=True) #change file name as needed
    print('Tscore done for dataset #' + str(datasetnum))
    print()
    return dfstat

def querygen(df, querysetsize, datasetnum):
    dfquery = pd.DataFrame()
    dfhighest = df['Tscore'].nlargest(querysetsize)
    dflowest = df['Tscore'].nsmallest(querysetsize)
    highestlist = pd.Index.tolist(dfhighest.index)
    lowestlist = pd.Index.tolist(dflowest.index)
    
    highestlist.insert(0, ('dfhighest' + str(datasetnum)))
    highestlist.insert(1, 'spacer')
    lowestlist.insert(0, ('dflowest' + str(datasetnum)))
    lowestlist.insert(1, 'spacer')
    
    highestSeries = pd.Series(highestlist)
    lowestSeries = pd.Series(lowestlist)
    dfquery = pd.concat([highestSeries, lowestSeries], axis = 1)
    dfquery = dfquery.T
    print('Query set generation done!')
    print()
    return dfquery

def prerankGSEA(refsigdf, datasetnum):
    outputdir = 'GSEAGeneID' + str(datasetnum)
    preresult = gp.prerank(refsigdf, 'Querysetdata.gmt', outdir=outputdir, pheno_pos='Pos', pheno_neg='Neg',
            min_size=15, max_size=500, permutation_num=1000, weighted_score_type=1,
            ascending=False, processes=1, figsize=(6.5,6), format='jpg',
            graph_num=20, no_plot=False, seed=None, verbose=False)
    preresult.run()

    terms = preresult.res2d.index
    i = 0
    for iter in terms:
        gseaplot(rank_metric=preresult.ranking, term=terms[i], **preresult.results[terms[i]])
        i = i + 1


df1 = pd.read_csv('GSE47960probes.txt', low_memory=False, delimiter = "\t") #change file name as needed
df2 = pd.read_csv('GSE47961probes.txt', low_memory=False, delimiter = "\t") #change file name as needed

normalizeddf = zscore(df1, 1)
controlcols = [67, 68, 69] #change column indices as needed for df1
experimentalcols = [137, 138, 139, 140] #change column indices as needed for df1
refsigdf1 = Tscore(normalizeddf, controlcols, experimentalcols, 1)

normalizeddf = zscore(df2, 2)
controlcols = [34, 35, 36] #change column indices as needed for df2
experimentalcols = [8, 9, 10, 11] #change column indices as needed for df2
refsigdf2 = Tscore(normalizeddf, controlcols, experimentalcols, 2)

dfquery1 = querygen(refsigdf1, 500, 1) #change 500 to needed query set size
dfquery2 = querygen(refsigdf2, 500, 2) #change 500 to needed query set size
querydf = pd.concat([dfquery1, dfquery2], axis = 0)
querydf.to_csv('Querysetdata.gmt', sep='\t', index=False)

prerankGSEA(refsigdf1, 1)
prerankGSEA(refsigdf2, 2)