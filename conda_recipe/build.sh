#!/bin/bash

#build.sh must be adapted to the logic in /src
#this small template is suitable for single file scripts
mkdir -p $PREFIX/bin/
cp get_fasta.py $PREFIX/bin/get_fasta
chmod +x $PREFIX/bin/
