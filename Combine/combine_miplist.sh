#!/bin/bash

#for filename in ~/Desktop/$1/Output/*miplist.txt*
for filename in ~/Desktop/$1/*miplist.txt*
do
    cp $filename ~/Desktop/Combine/
done

output=""
for filename in ~/Desktop/Combine/*miplist.txt*
do
	output="$output $filename"
done

python combine_miplist.py $output

sort -nk2,2 ~/Desktop/Combine/miplist.combined.txt | sort -nk1,1 > ~/Desktop/Combine/$1.miplist.combined.txt
	
for filename in ~/Desktop/Combine/*miplist.txt*
do
rm $filename
done