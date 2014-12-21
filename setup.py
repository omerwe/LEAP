"""
file to set up python package, see http://docs.python.org/2/distutils/setupscript.html for details.
"""


import platform
import os
import sys
import shutil
import setuptools

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.clean import clean as Clean


try: import numpy	
except Exception:
	print "numpy needed for installation, please install numpy first"
	sys.exit()
try: import scipy
except Exception:
	print "scipy needed for installation, please install numpy first"
	sys.exit()

def readme():
    with open('leap/README.md') as f:
       return f.read()

# set up macro
if "win" in platform.system().lower():
    macros = [("_WIN32", "1")]
else:
    macros = [("_UNIX", "1")]


#python setup.py sdist bdist_wininst upload
setup(
    name='leap_gwas',
    version='0.1.4.1',
    description='Liability Estimation in Case Control Studies',
    long_description=readme(),
    keywords='gwas bioinformatics LMMs MLMs',
    url="https://github.com/omerwe/LEAP",
    author='Omer Weissbrod',
    author_email='omerw@cs.technion.ac.il',
    license='Apache 2.0',    
    packages=["leap", "leap/leap", "leap/dataset1", "leap/results_gold"],
    include_package_data=True,
    install_requires = ['numpy', 'scipy', 'fastlmm', 'sklearn'],
    #zip_safe=False,
    # extensions        
	include_dirs = [numpy.get_include()]
  )	   
	
