#/.../.../mipgen_analysis.sh folder read1.fq read2.fq prefix

threads=4
script_dir=~/MIPGEN/tools/
#barcodes=~/MIPGEN/reference/barcode_blank.txt
barcodes=~/Desktop/$1/barcode_file.txt
mips=~/MIPGEN/reference/ordered_polyg_8N.design.subset.txt
prefix=$4
pear_read1=~/Desktop/$1/$2
pear_read2=~/Desktop/$1/$3
ext_tag_size=4
lig_tag_size=4
fqprefix=~/Desktop/$1/${prefix}
bamprefix=${fqprefix}.indexed.sort
cut_read1=${fqprefix}.assembled.fastq
#gref=~/MIPGEN/reference/human_g1k_v37.fasta
gref=~/MIPGEN/reference/synthetic_reference.fasta

pear -j $threads -f $pear_read1 -r $pear_read2 -o $fqprefix &&
python ${script_dir}mipgen_fq_cutter_se.py $cut_read1 -tb $barcodes -m ${lig_tag_size},${ext_tag_size} -o $fqprefix &&
~/tools/bin/bwa-0.7.12/bwa mem -t $threads $gref ${fqprefix}.indexed.fq > ${fqprefix}.indexed.sam &&
samtools view -bS ${fqprefix}.indexed.sam | samtools sort - $bamprefix &&
samtools view -h $bamprefix.bam | python ${script_dir}mipgen_smmip_collapser.py 8 $bamprefix.collapse -m $mips -f 1 -c -r -b $barcodes -w -s &&
echo "analysis commands have terminated (successfully or otherwise)"

ulimit -n 10000
samtools view -S ${fqprefix}.indexed.sort.collapse.all_reads.uncollapsed.sam | python ${MipDir}explain_cigars.py $barcodes ${MipDir}subset_mip_scanstart_relstart_len_type.8N.txt ${fqprefix}.explained 1
sort -nk2,2 ${fqprefix}.explained.shortform.txt | sort -nk1,1 > ${fqprefix}.explained.shortform.sorted.txt
#echo 'mip       tract_length        count       fraction' | cat - ${fqprefix}.explained.shortform.sorted.txt  > temp && mv temp ${fqprefix}.explained.shortform.sorted.txt

#mkdir ~/Desktop/$1/$prefix
#for file in ~/Desktop/$1/*explained*
#do
#mv $file ~/Desktop/$1/$prefix
#done

#Rscript ${MipDir}MipRscript.r ${prefix}
echo "Done!"