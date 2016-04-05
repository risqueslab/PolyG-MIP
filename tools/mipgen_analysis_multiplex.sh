#/.../.../mipgen_analysis_multiplex.sh folder read1.fq read2.fq readi.fq barcode prefix

set -e
set -o pipefail

threads=4
script_dir=~/MIPGEN/tools/
indexread=~/Desktop/$1/$4
barcodes=~/Desktop/$1/$5
mips=~/MIPGEN/reference/ordered_polyg_8N.design.subset.txt
prefix=$6
pe_premerge_read1=~/Desktop/$1/$2
pe_premerge_read2=~/Desktop/$1/$3
pear_read1=~/Desktop/$1/${prefix}.r1.indexed.fq
pear_read2=~/Desktop/$1/${prefix}.r2.indexed.fq
ext_tag_size=4
lig_tag_size=4
fqprefix=~/Desktop/$1/${prefix}
bamprefix=${fqprefix}.indexed.sort
cut_read1=${fqprefix}.assembled.fastq
gref=~/MIPGEN/reference/human_g1k_v37.fasta
#gref=~/MIPGEN/reference/synthetic_reference.fasta

python ${script_dir}mipgen_fq_cutter_pe.py $pe_premerge_read1 $pe_premerge_read2 -tb $barcodes -i $indexread -j 10 -o $fqprefix &&
pear -j $threads -f $pear_read1 -r $pear_read2 -o $fqprefix &&
python ${script_dir}mipgen_fq_cutter_se.py $cut_read1 -tb $barcodes -m ${lig_tag_size},${ext_tag_size} -o $fqprefix &&
~/tools/bin/bwa-0.7.12/bwa mem -t $threads $gref ${fqprefix}.indexed.fq > ${fqprefix}.indexed.sam &&
samtools view -bS ${fqprefix}.indexed.sam | samtools sort - $bamprefix &&
samtools view -h $bamprefix.bam | python ${script_dir}mipgen_smmip_collapser.py 8 $bamprefix.collapse -m $mips -f 1 -c -r -b $barcodes -w -s &&
echo "analysis commands have terminated (successfully or otherwise)"

samtools view -S ${fqprefix}.indexed.sort.collapse.all_reads.uncollapsed.sam | python ${MipDir}explain_cigars.py $barcodes ${MipDir}subset_mip_scanstart_relstart_len_type.8N.txt ${fqprefix}.explained 1
sort -nk3,3 ${fqprefix}.explained.shortform.txt | sort -nk1,1 |sort -sk2,2 > ${fqprefix}.explained.shortform.sorted.txt
#echo 'mip  sample  tract_length        count       fraction' | cat - ${fqprefix}.explained.shortform.sorted.txt  > temp && mv temp ${fqprefix}.explained.shortform.sorted.txt

#mkdir ~/Desktop/$1/$prefix
#for file in ~/Desktop/$1/*explained*
#do
#mv $file ~/Desktop/$1/$prefix
#done

#Rscript ${MipDir}MipRscript.r ${prefix}
echo "Done!"

