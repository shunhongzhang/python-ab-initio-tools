#!/usr/bin/python

import os
import numpy as np

def parse_pwi(filename='scf.in'):
    ibrav = os.popen('grep ibrav '+filename).read()
    ibrav = int(ibrav.split('=')[-1].rstrip('\n').rstrip(','))
    print 'ibrav= {0:1d}'.format(ibrav)
    celldm=np.zeros(7,float)
    for i in range(1,7):
        get_celldm='grep "celldm('+str(i)+')" '+filename
        try:
           celldm[i] = float(os.popen(get_celldm).read().split('=')[-1].rstrip('\n').rstrip(','))
           print celldm[i]
        except:
           pass
    if ibrav==0:
       print 'free lattice read'
    elif:
       ibrav==1 or ibrav==

parse_pwi()
