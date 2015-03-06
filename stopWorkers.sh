#!/bin/bash
echo killing all current workers in current directory
sudo kill `ps aux | grep $(pwd) | grep -v grep | tr -s ' ' | cut -f2 -d ' '`
