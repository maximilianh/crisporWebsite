#!/bin/bash
# Input file list
INPUT_FILES=(logdbl.cpp)
# Output file prefix
PREFIX=logdbl
# STDIN input when running the program(s)
STDIN_DATA="$1"
# Custom g++ #defines -- One array element for each program.
# e.g. BUILD_FLAGS=("-DFRED -DALICE"  "-w -DBOB"  "-std=c++11 -D BANANA")
# Not-yet-implimented: BUILD_FLAGS='-DD#' is A shortcut for using -DD{N} (where N is the 1-based program number)
BUILD_FLAGS=("${@:2}");  BUILD_FLAGS="${BUILD_FLAGS[@]/#/-DD}"

# Build operations R=Run  A=diff-asm  S=diff-src
OPS=RAS
# Common g++ Flags for all programs
CPP_FLAGS=-O3
# flags for objdump
OBJ_FLAGS='-drwCS -Mintel'

source run-test-builds.sh
