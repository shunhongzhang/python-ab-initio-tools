#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       setpu.py                                                     #
#  Dependence: the whole tool-kit                                           #
#  Usage:      install the files as a lib and generate excutables           #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Apr 15, 2017                                                 #
#                                                                           #
#===========================================================================#


from distutils.core import setup
#from setuptools import setup

my_modules=['crystal_structure','utils','parse','write','match_latt','v2qe','v2openmx','dos','kpdos','bands','arguments','pythtb_plus','grid']


setup(name='v2openmx',
      version='1.0.0',
      author='Shunhong Zhang',
      author_email='zhangshunhong@tsinghua.edu.cn',
      url='https://to_be_posted',
      download_url='https://on_request',
      keywords='pre- and post-processing for ab initio compuational materials science',
      py_modules=my_modules,
      license="gpl-3.0",
      description="pre- and post-processing for ab initio compuational materials science",
      long_description="A collection of scripts for pre- and post-processing of files related to ab initio DFT calculations, interfacing with codes like VASP, Quantum ESPRESSO, OPENMX, and phonopy",
      platforms=['LINUX'],
      #install_requires=['numpy','matplotlib','yaml','termcolor'],
      )

import os
pwd=os.popen('pwd').read().rstrip('\n')
if not os.popen('grep "'+pwd+'" ~/.bashrc|grep PATH').read():
   os.system('echo "export PATH=\$PATH:\c"|cat >>~/.bashrc')
   os.system('pwd|cat >>~/.bashrc')
   os.system('echo "export PATH=\$PATH:"'+pwd+'/commands_abbr|cat >>~/.bashrc')
   os.system('echo "export PYTHONPATH=\$PYTHONPATH:"'+pwd+'/lib/python|cat >>~/.bashrc')


print '\nsuccessfully installed, now test\n'


import importlib

for mod in my_modules:
    try:
       importlib.import_module(mod)
       print 'import module {0:20s} successfully'.format(mod)
    except:
       print 'import module {0:20s} failed!'.format(mod)

print '\nplease report to the author if you failed in importing some modules.'
print 'thank you for your feedback!'
