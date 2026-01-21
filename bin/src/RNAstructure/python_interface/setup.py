#!/usr/bin/env python

"""
setup.py file for RNAstructure python interface
"""

import os
import platform
from distutils.core import setup, Extension
from distutils.sysconfig import get_config_vars

# This removes the annoying -Wstrict_prototypes Warning
(opt,) = get_config_vars('OPT')
os.environ['OPT'] = " ".join(flag for flag in opt.split() if flag != '-Wstrict-prototypes')

# Read rna_sources.h and include each line as a source file name.
# (Ignore blank lines and lines starting with '#')
with open('rna_sources.h') as f:
  sources = [line.rstrip('\n') for line in f if line != '' and line[0] != '#' ]

CXXFLAGS = ['-O3', '-w']
LDFLAGS = [ '-lpthread' ]

# static_args are only used on Windows.
#static_args = ['-static', '-static-libgcc', '-static-libstdc++']

#if os.name == 'nt' or "CYGWIN" in platform.system():
#  LDFLAGS.extend(static_args)

RNAstructure_module = Extension('_RNAstructure_wrap',
								sources=sources,
								extra_compile_args=CXXFLAGS,
                extra_link_args = LDFLAGS,
								)

# Change PWD to RNAstructure root directory, so 
# relative paths to cpp files will work.
script_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(os.path.dirname(script_dir))

setup (name = 'RNAstructure_wrap',
       version = '1.0',
       author      = "Michael Sloma",
       description = """Python interface for the RNA library""",
       ext_modules = [RNAstructure_module],
       py_modules = ["RNAstructure_wrap", "Error_handling"],
       )
