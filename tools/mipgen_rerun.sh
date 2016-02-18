#Alexander M. West
#5/29/2015
#generate new longform and shortforms with already processed sam files
#syntax: sh mipgen_rerun.sh folder

set -e
set -o pipefail
ulimit -n 10000

for filename in ~/Desktop/$1/*.indexed.sort.collapse.all_reads.uncollapsed.sam

do
samtools view -S ~/Desktop/$1/$(basename ${filename/.*/}).indexed.sort.collapse.all_reads.uncollapsed.sam | python ~/MIPGEN/tools/explain_cigars.py ~/Desktop/$1/barcode_file.txt ~/MIPGEN/tools/subset_mip_scanstart_relstart_len_type.8N.txt ~/Desktop/$1/$(basename ${filename/.*/}).explained 1
sort -nk2,2 ~/Desktop/$1/$(basename ${filename/.*/}).explained.shortform.txt | sort -nk1,1 > ~/Desktop/$1/$(basename ${filename/.*/}).explained.shortform.sorted.txt
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
for file in ~/Desktop/$1/*_out*
do
mv $file ~/Desktop/$1/Output
done