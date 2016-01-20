#!/bin/bash
echo killing all current workers in current directory
dirName=`pwd`
dirName=`basename $dirName`
echo directory: $dirName
kill `ps aux | grep $dirName | grep -v grep | tr -s ' ' | cut -f2 -d ' '`
