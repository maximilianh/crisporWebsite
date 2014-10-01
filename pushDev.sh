#!/bin/bash
mv crispor.cgi crispor.cgi.old
cp crisporDev.cgi crispor.cgi

mv bin bin.old
cp -R binDev/ bin/
