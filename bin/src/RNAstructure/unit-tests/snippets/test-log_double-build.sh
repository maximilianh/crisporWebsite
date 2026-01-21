#!/bin/bash
# Input file list
INPUT_FILES=(test-log_double.cpp)
# Output file prefix
PREFIX=tld_${3:-2}
# STDIN input when running the program(s)
STDIN_DATA='' #"100 11.25 1E6" #"$2"
# Custom g++ #defines -- One array element for each program.
# e.g. BUILD_FLAGS=("-DFRED -DALICE"  "-w -DBOB"  "-std=c++11 -D BANANA")
# Not-yet-implimented: BUILD_FLAGS='-DD#' is A shortcut for using -DD{N} (where N is the 1-based program number)
BUILD_FLAGS=(-DD1 -DD2) #("${@:3}");  BUILD_FLAGS="${BUILD_FLAGS[@]/#/-DD}"

DEF_OPS=BARS
DEF_OBJ=-drwCS

# Build operations R=Run  A=diff-asm  S=diff-src
OPS="${1:-$DEF_OPS}" # Default is to perform all operations: 'BASR'
# Common g++ Flags for all programs
CPP_FLAGS="-O3 -pedantic -DCOUNT=${3:--1}"
# flags for objdump
OBJ_FLAGS="${2:-$DEF_OBJ}"  # Default if empty: '-drwCS', See `man objdump` for details.


source run-test-builds.sh
