#!/bin/bash

for filename in ~/Desktop/$1/*notes.txt*
do
    cp $filename ~/Desktop/Combine/
done

output=""
for filename in ~/Desktop/Combine/*notes.txt*
do
	output="$output $filename"
done

python combine_notes.py $output

#sort -nk2,2 ~/Desktop/Combine/notes.combined.txt | sort -nk1,1 > ~/Desktop/Combine/$1.notes.combined.txt

for filename in ~/Desktop/Combine/*notes.txt*
do
rm $filename
done

	