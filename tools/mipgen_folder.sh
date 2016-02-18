#!/bin/bash
ulimit -n 10000
for filename in ~/Desktop/$1/*R1*
do
    sh mipgen_analysis.sh $1 $(basename $filename) $(basename ${filename/R1/R2}) $(basename ${filename/R1*/Analysis})
done

#script to move outputs to their own folder
mkdir ~/Desktop/$1/Output
for file in ~/Desktop/$1/*longform*
do
mv $file ~/Desktop/$1/Output
done
for file in ~/Desktop/$1/*sorted*
do
mv $file ~/Desktop/$1/Output
done
for file in ~/Desktop/$1/*histogram*
do
mv $file ~/Desktop/$1/Output
done
for file in ~/Desktop/$1/*miplist*
do
mv $file ~/Desktop/$1/Output
done
