#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
from twoLay import twolay
import pickle
import random
import itertools
import time
import pickle
import sys
from scipy.interpolate import griddata


import dask
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster

if __name__ == "__main__":

    cluster = LocalCluster()
    client = Client(cluster)
    
    
    species = ['OAH1','SAN1', 'AUS1','ICE1','KUW1','NIG1','SAH1','DET']
    nprimepath = 'stramski_2007_mineral_nprime.csv'
    nprime = pd.read_csv(nprimepath,index_col=0)
    
    # phytos = phytodata.Species
    species_dict = {}
    
    count = 0
    name = 0
    for p in species:
        if count%1 == 0:
            name = name+1
            species_dict[name] = []
            species_dict[name].append(p)
        else:
            species_dict[name].append(p)
        count+=1

    for key in species_dict.keys():
        print('On node {} are processed: {}'.format(key, species_dict[key]))
    key = int(sys.argv[1])
    print('the key is: {}'.format(key))

    parameters = ['Qc',
                  'Sigma_c',
                  'cstar',
                  'Qb',
                  'Sigma_b',
                  'bstar',
                  'Qa',
                  'Sigma_a',
                  'astar',
                  'Qbb',
                  'Sigma_bb',
                  'bbstar',
                  'VSF',]
    
    # wavelength range and resolution 
    #(changing this changes your interp value when normalising kshell)
    l = np.arange(.4, .9025, .0025).astype(np.float32) 
    
    
    count = 0
    data = {}
    for species in species_dict[key]:
        #start timer
        start = time.time()

        #add phyto to dictionary
        data[species] = {}

        
        if species == 'DET':
            nreal = [1.03, 1.05]
            jexp = [3.4, 4, 4.6]
            dmax = [10.1, 50.1, 100.1]
            rho = [.3e6, .5e6, .7e6]  
            kcore = 0.010658 * np.exp(-0.007186* (l*1000)) #Stramski 2001 
            kcore = kcore.squeeze()
        else:    
            nreal = [1.1, 1.4]
            jexp = [3.4, 4, 4.6]
            dmax = [10.1, 50.1, 100.1]
            rho = [.75e6, 2e6, 4e6]
            kcore = nprime[species].values.copy()  
            im_wv = nprime.index.values / 1000
            last = kcore[-1:]
            kcore = griddata(im_wv, kcore, l, 'linear',fill_value=last)
            kcore = kcore.squeeze()
        
        #create iterlist
        iterlist = []
        for it in itertools.product(nreal, jexp, dmax, rho):
            run = {}
            run['nreal'] = it[0]
            run['jexp'] = it[1]
            run['dmax'] = it[2]
            run['rho'] = it[3]
            iterlist.append(run)

        #create dictionary entries
        for i in range(len(iterlist)):
            rname = '{}_{:.2f}_{:.2f}_{:.2f}_{}'.format(species,
                                                        iterlist[i]['nreal'],
                                                        iterlist[i]['jexp'],
                                                        iterlist[i]['dmax'],                                                                  
                                                        iterlist[i]['rho'],
                                                        )   
            data[species][rname] = iterlist[i]        
    
    
        for rname in data[species].keys():
            print(rname)

            # RUN twoLay
            result0 = dask.delayed(twolay)(l, 
                                           kcore, 
                                           data[species][rname]['nreal'], 
                                           data[species][rname]['jexp'], 
                                           data[species][rname]['dmax'], 
                                           data[species][rname]['rho'],)
            for p in result0.keys():
                data[species][rname][p] = result0[p]
            
        data[species] = dask.compute(data[species])[0]
        
        #end timer
        end = time.time()
        tleng = end - start
        data[species]['time'] = tleng    
    
    for species in data.keys():
        print('Species {} ran for {:.2f} seconds.'.format(species, data[species]['time']))

    with open('napdata_pbs_01.dat', 'wb') as picklefile:
        pickle.dump(data, picklefile)    
    
    
    