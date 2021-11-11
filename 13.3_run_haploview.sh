#!/bin/sh

for i in `ls *.haps`; do
    fullname=$(basename "$i");
    filename="${fullname%.*}";
    mkdir -p res$filename;
    java -jar ./Haploview.jar -n -memory 200000 -log $filename.log -haps $i -info $filename.info -dprime -out $filename.res -minGeno 0.2 -missingCutoff 0.7 -hwcutoff 0 -blockoutput GAB; #in Kb
    mv $filename.res.LD $filename.res.GABRIELblocks res$filename;
done
