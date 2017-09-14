#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       readphonon.py                                                #
#  Dependence: none                                                         #
#  Usage:      parse data from band.yaml and plot the phonon spectra        #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Apr 22, 2017                                                 #
#                                                                           #
#===========================================================================#

import numpy as np
try:
    import yaml
except ImportError:
    sys.exit("You need to install python-yaml")
try:
    from yaml import CLoader as Loader
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

def parse_band_yaml(filename):
    data = yaml.load(open(filename), Loader=Loader)
    frequencies = []
    distances = []
    labels = []
    for j, v in enumerate(data['phonon']):
        if 'label' in v:
            labels.append(v['label'])
        else:
            labels.append(None)
        frequencies.append([f['frequency'] for f in v['band']])
        distances.append(v['distance'])
    return (np.array(distances), np.array(frequencies),
            data['segment_nqpoint'],labels)

distance,freq,segment_nqpoint,labels=parse_band_yaml('band.yaml')
nqpt,nband=freq.shape
print 'distance of high symmetry q-points:'
for iseg in range(len(segment_nqpoint)):
    print labels[iseg],distance[np.cumsum(segment_nqpoint)[iseg]-1]


fw = open('freq.dat','w')
for iband in range(nband):
    print>>fw,'\n'.join('{0:10.7f} {1:10.7f}'.format(dist,fre) for dist,fre in zip(distance,freq[:,iband])),'\n'
fw.close()
fw=open('freq-sort.dat','w')
for iqpt in range(nqpt):
    print>>fw,distance[iqpt],
    for iband in range(nband):
        print>>fw,'{0:10f}'.format(freq[iqpt,iband]),
    print>>fw,''
fw.close()
