#!/usr/bin/python

# This script can be used to show the value of a python distutils config variable.
# Author: Richard M. Watson (2018-2019)

# Usage:
# get-config.py [list|all|none|VARIABLES...]
#    list -- Print a list of all config variable names (there are a lot!)
#    all  -- Print all variable names with their values.
#    VARIABLES -- The name(s) of one or more variables, separated by spaces.
#      Variable names are usually upper-case.
#      If more than one variable name is listed, their values are printed on
#      separate lines.
#
#    If a variable name is invalid, nothing is printed, and the exit code will be 1.
#
# Sometimes variables will be valid for one version of python and not another.
# So there is a special format used to get the value of the first VALID name
# in a list. Separate multiple names with "/" (a slash).
# for example:  `get-config.py LDVERSION/VERSION` -- This would print the value
#   of LDVERSION if it is a valid name, or VERSION otherwise.

# Examples:
#   get-config.py list       -- Print a list of all config variable names (there are a lot!)
#   get-config.py all        -- Print all variable names with their values.
#   get-config.py VERSION    -- show the python version number e.g. 3.6
#   get-config.py EXT_SUFFIX -- show the file suffix for python extensions on this OS. (python 3 only)

import os
import distutils
import distutils.sysconfig
#from setuptools import setup
import sys
import datetime

# print the value of the variable name
def print_val(name):
  if name == 'all':
    print_all(True)
  elif name == 'list':
    print_all(False)
  elif '/' in name:
    return print_first_valid(name.split('/'))
  else:
    val=distutils.sysconfig.get_config_var(name)
    if val is None:
      return False
    print(val)
  return True # always return true, unless a branch above already returned.

# Print all variables (and their values, if showVals is true)
def print_all(showVals):
  vars=distutils.sysconfig.get_config_vars()
  for n in vars:
    if showVals:
      print(str(n) + ":\t" + str(vars[n]))
    else:
      print(str(n))
  return True

# Print the first valid name in the list
def print_first_valid(names):
  for name in names:
    val=distutils.sysconfig.get_config_var(name)
    if not val is None:
      print(val)
      return True
  return False


names = sys.argv[1:] # skip first two args (python and script path)
if len(names)==0:
  print("ERROR -- at least one argument is required -- a distutils sysconfig variable name (e.g. VERSION).")
  quit(1)

retval=0
for name in names:
  if not print_val(name):
    retval = 1

quit(retval)
