#!/bin/bash

#Useage: combine_shortform.sh <target folder>

for filename in ~/Desktop/$1/*sorted*
#for filename in ~/Desktop/$1/Output/*sorted*
do
cp $filename ~/Desktop/Combine/
done

output=""
for filename in ~/Desktop/Combine/*sorted*
do
	output="$output $filename"
done

python combine_shortform.py $output
#python combine_shortform_likemultiplex.py $output

sort -nk2,2 ~/Desktop/Combine/shortform.combined.txt | sort -nk1,1 > ~/Desktop/Combine/$1.shortform.combined.txt
	
for filename in ~/Desktop/Combine/*sorted*
do
rm $filename
done