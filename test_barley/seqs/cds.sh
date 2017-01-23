#!/bin/bash
for f in *.gz *.bz2
do
  echo "Processing $f file..."
  $GETHOMS/transcripts2cds.pl -X -n 22 $f &> log.$f
done
