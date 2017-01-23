#!/bin/bash
for f in *.gz *.bz2
do
  echo "Processing $f file..."
  $GETHOMS/transcripts2cdsCPP.pl -X -n 20 $f &> log.$f
done
