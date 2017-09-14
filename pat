#!/usr/bin/python

import numpy as np
import sys
import os
import argparse
import arguments

desc_str = '''python ab-initio tools (PAT)'''
parser = argparse.ArgumentParser(prog='pat', description = desc_str)
parser.add_argument('task',type=str,default=None,help='The task to execute')
parser.add_argument('input',type=str,default='pat.in',help='filename of the input')
arguments.add_io_arguments(parser)
arguments.add_fig_arguments(parser)
arguments.add_plot_arguments(parser)
args = parser.parse_args()

def parse_input(args):
    filinput=args.input
    f=open(filinput)
    data = np.loadtxt(f,comments='#',dtype=str)
    arguments = {}
    for idata in data:
        key   = idata.split('=')[0].lower()
        value = idata.split('=')[1]
        arguments.setdefault(key,value)
    return arguments


arguments=parse_input(args)
cmd=args.task
for tag in arguments.keys():
    cmd+=' --'+tag+'='+arguments[tag]
print 'executing the following command:\n',cmd
os.system(cmd)
