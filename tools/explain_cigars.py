# Maintained by Alexander M. West
# Version: Feb. 29, 2016

import sys
import os
sys.path.append("~/MIPGEN/tools")
from genome_sam_collapser import *
from consensusMaker import *
import re
from collections import defaultdict
counts = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : 0)))
totals = defaultdict(lambda : defaultdict(lambda : 0))
topdepth = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : (0,0))))
alldepth = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : 0))))
barcode_lookup = {}
barcode_labels = {}
barcode_file = open(sys.argv[1])
mip_scanstart_relstart_len = open(sys.argv[2])
prefix = sys.argv[3]
count_threshold = int(sys.argv[4])
read_agreement_percentage_threshold = float(0)/10 #percentage of reads in a tag that must agree to call the tag
minimum_read_threshold = 2
minimum_tag_threshold = 1

allele_out = open(prefix + ".allele_out.txt", 'w')
cigar_out = open(prefix + ".cigar_out.txt", 'w')
longform = open(prefix + ".longform.txt", 'w')
shortform = open(prefix + ".shortform.txt", 'w')
debug_cigars = open(prefix + ".debug_cigars.txt", 'w')
histogram = {}
miplist = {}
#mip_in_process = {}
barcode_tag_list = {}

mip_key_first_coordinate_lookup = {}
mip_to_info = {}
mip_sample_alleles = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : 0)))

#initialize barcode information
for line in barcode_file:
  label, seq = line.rstrip().split()
  barcode_labels[seq] = label
  bases = "ATCGN"
  for i in range(len(seq)):
    native_seq = list(seq)
    for j in range(5):
      mutated_seq = native_seq
      mutated_seq[i] = bases[j]
      barcode_labels["".join(mutated_seq)] = barcode_labels[seq]
  histogram[label] = open(prefix + "." + label + ".histogram.txt", 'w')
  histogram[label].write("Reads/Tag\tFrequency\n")
  miplist[label] = open(prefix + "." + label + ".miplist.txt", 'w')
  miplist[label].write("MIP\tReads\tTags < 3\tTags > 2\n") #Unused Tags\tUsed Tags\n")  Unused labels for columns 1 and 2.
barcode_file.close()

#initialize mip location information
for line in mip_scanstart_relstart_len:
  mip_key, scanstart, relstart, tractlen, type = line.rstrip().split()
  mip_key_first_coordinate_lookup[scanstart] = mip_key
  mip_to_info[mip_key] = (int(relstart), int(tractlen))
mip_scanstart_relstart_len.close()
mip_key_second_coordinate_lookup = {}
cigar_pattern = re.compile("([\d]+)([A-Z])")

#iterate through reads, sort into tags and inventory
for sam_line in sys.stdin:
  current_read = sam_read(sam_line, mip_key_first_coordinate_lookup, mip_key_second_coordinate_lookup)
  debug_cigars.write(current_read.mip_key + "\t" + current_read.mtag + "\t" + current_read.seq + "\t")
  if current_read.barcode not in barcode_labels:
    continue
  relstart, tractlen = mip_to_info[current_read.mip_key]
  debug_cigars.write(str(relstart) + "\t" + str(tractlen) + "\t")
  parsed_cigar = []
  for pair in cigar_pattern.findall(current_read.cigar):
    parsed_cigar.append([int(pair[0]), pair[1]])
  debug_cigars.write(str(parsed_cigar) + "\t")
  current_read.trim_scan_sequences(relstart - 1, relstart + tractlen -2, parsed_cigar)
  debug_cigars.write(current_read.seq + "\n")
		# -2 to correct for awk index and get adjacent base, +1 to get adjacent base,
		# -1 to get to last scan, -1 to correct for awk index
		#  print current_read.mip_key, current_read.seq, current_read.cigar
        # +1, -1 to remove all adjacent bases
  #sequence_path = prefix.rpartition("/")[0] + "/Sequences/" + barcode_labels[current_read.barcode] + "/"
  sequence_path = prefix.rpartition("/")[0] + "/Sequences/" + prefix.rpartition("/")[2].rpartition("_")[0] + "/"
  if not os.path.exists(sequence_path):
    os.makedirs(sequence_path)
  file = open(sequence_path + current_read.mip_key.split("/")[0] + "." + current_read.mtag + ".txt", 'a')
  file.write("> " + current_read.mip_key + "\t" + current_read.mtag + "\t" + current_read.cigar + "\n" + current_read.seq + "\n")
  file.close()

  obs_tract = current_read.cigar_seq_length()
  refoverlap = current_read.cigar_length()
  if len(current_read.seq) == 0:
    continue
  #line below used to be tractlen + n, where n is the number of bases around the target kept after trim_scan_sequences
  if refoverlap == tractlen + 0:
    if re.search("(G{9}G*|C{9}C*)", current_read.seq):
      mip_sample_alleles[current_read.mip_key][current_read.barcode][current_read.seq + " homo"] += 1
    else:
      mip_sample_alleles[current_read.mip_key][current_read.barcode][current_read.seq] += 1
    current_depth = alldepth[current_read.mip_key][barcode_labels[current_read.barcode]][current_read.mtag][obs_tract] + 1
    alldepth[current_read.mip_key][barcode_labels[current_read.barcode]][current_read.mtag][obs_tract] = current_depth
    if current_depth >= count_threshold:
      if topdepth[current_read.mip_key][current_read.barcode][current_read.mtag] == (0, 0):
        counts[current_read.mip_key][barcode_labels[current_read.barcode]][obs_tract] += 1
        totals[current_read.mip_key][barcode_labels[current_read.barcode]] += 1
        topdepth[current_read.mip_key][current_read.barcode][current_read.mtag] = (obs_tract, 1)
      elif topdepth[current_read.mip_key][current_read.barcode][current_read.mtag][1] <= current_depth:
        old_tract = topdepth[current_read.mip_key][current_read.barcode][current_read.mtag][0]
        old_depth = topdepth[current_read.mip_key][current_read.barcode][current_read.mtag][1]
        if old_tract != obs_tract:
          counts[current_read.mip_key][barcode_labels[current_read.barcode]][old_tract] -= 1
          counts[current_read.mip_key][barcode_labels[current_read.barcode]][obs_tract] += 1
        topdepth[current_read.mip_key][current_read.barcode][current_read.mtag] = (obs_tract, old_depth+1)

mip_allele_count = {}   #a dictionary of alleles and frequency of tags
allele_tag_total = 0    #the number of tags for a given allele?
tag_read_tally = {}
mip_read_tally = {}

#adding headers to longform file
#DN: longform.write("MIP" + "\t" + "SAMPLE" + "\t" + "TAG" + "\t" + "SEQUENCE" + "\t" + "READS" + "\n")

#iterate through all reads and perform appropriate methods
for mip, bar_dict in alldepth.iteritems():
  for bar, mtag_dict in bar_dict.iteritems():
    for mtag, tract_dict in mtag_dict.iteritems():
      total_tcount = 0;
      for tract, tcount in tract_dict.iteritems():
#allele count level
        #print tract, tcount
        total_tcount += tcount
        if (mip, bar) in mip_read_tally.keys():
            mip_read_tally[(mip, bar)] = (mip_read_tally[(mip, bar)][0] + tcount, mip_read_tally[(mip, bar)][1], mip_read_tally[(mip, bar)][2], mip_read_tally[(mip, bar)][3], mip_read_tally[(mip, bar)][4])
        else:
            mip_read_tally[(mip, bar)] = (tcount, 0, 0, 0, 0)
#molecular tag level
      #print mtag_dict
      #print mtag, tract_dict
      #print total_tcount
      if total_tcount > 2:
          mip_read_tally[(mip, bar)] = (mip_read_tally[(mip, bar)][0], mip_read_tally[(mip, bar)][1], mip_read_tally[(mip, bar)][2], mip_read_tally[(mip, bar)][3], mip_read_tally[(mip, bar)][4] + 1)
      else:
          mip_read_tally[(mip, bar)] = (mip_read_tally[(mip, bar)][0], mip_read_tally[(mip, bar)][1], mip_read_tally[(mip, bar)][2], mip_read_tally[(mip, bar)][3] + 1, mip_read_tally[(mip, bar)][4])
      if (total_tcount, bar) in tag_read_tally.keys():
        tag_read_tally[(total_tcount, bar)] = tag_read_tally[(total_tcount, bar)] + 1
      else:
        tag_read_tally[(total_tcount, bar)] = 1
      if total_tcount > 2:
          #sequence_path = prefix.rpartition("/")[0] + "/Sequences/" + bar + "/"
        sequence_path = prefix.rpartition("/")[0] + "/Sequences/" + prefix.rpartition("/")[2].rpartition("_")[0] + "/"
        os.system("linsi --allowshift --kappa 1 --quiet --thread 4 " + sequence_path + mip.split("/")[0] + "." + mtag + ".txt > " + sequence_path + mip.split("/")[0] + "." + mtag + ".aligned.txt")
        sequence, sequence_info = consensusMaker(sequence_path + mip.split("/")[0] + "." + mtag + ".aligned.txt")
        sequence = sequence.replace("-", "")
        longform.write(mip + "\t" + bar + "\t" + mtag + "\t" + sequence + "\t" + sequence_info + "\n")
        if not sequence == "":
          if sequence in mip_allele_count.keys():
            mip_allele_count[sequence] = mip_allele_count[sequence] + 1
          else:
            mip_allele_count[sequence] = 1
#barcode level
    for allele in mip_allele_count:
      if mip_allele_count[allele] > minimum_tag_threshold:
        allele_tag_total += mip_allele_count[allele]
        mip_read_tally[(mip, bar)] = (mip_read_tally[(mip, bar)][0], mip_read_tally[(mip, bar)][1], mip_read_tally[(mip, bar)][2] + mip_allele_count[allele], mip_read_tally[(mip, bar)][3], mip_read_tally[(mip, bar)][4])
      else:
        mip_read_tally[(mip, bar)] = (mip_read_tally[(mip, bar)][0], mip_read_tally[(mip, bar)][1] + mip_allele_count[allele], mip_read_tally[(mip, bar)][2], mip_read_tally[(mip, bar)][3], mip_read_tally[(mip, bar)][4])
    for allele in mip_allele_count:
      if mip_allele_count[(allele)] > minimum_tag_threshold:
          #shortform.write(mip + "\t" + bar + "\t" + str(allele[0]) + "\t" + str(allele[1]) +"\t" + str(mip_allele_count[allele]) + "\t" + str(mip_allele_count[allele] / float(allele_tag_total)) + "\n")
          shortform.write(mip + "\t" + bar + "\t" + str(allele) + "\t" + str(len(allele)) + "\t" + str(mip_allele_count[allele]) + "\t" + str(mip_allele_count[allele] / float(allele_tag_total)) + "\n")
    mip_allele_count = {}
    allele_tag_total = 0
#mip level output

#output the dictionary of family sizes to experiment specific text files
for value in sorted(tag_read_tally, key=lambda read: read[0]):
    histogram[value[1]].write(str(value[0]) + "\t" + str(tag_read_tally[value]) + "\n")

#currently not printing value 1 or 2.
for value in sorted(mip_read_tally, key=lambda read: read[0].partition(":")):
    miplist[value[1]].write(value[0] + "\t" + str(mip_read_tally[value][0]) + "\t" + str(mip_read_tally[value][3])+ "\t" + str(mip_read_tally[value][4]) + "\n")

#close all
longform.close()
shortform.close()
debug_cigars.close()
for label in histogram:
    histogram[label].close()
for label in miplist:
    miplist[label].close()
for mip, sample_dict in mip_sample_alleles.iteritems():
  for sample, allele_dict in sample_dict.iteritems():
    for allele, c in allele_dict.iteritems():
      allele_out.write(mip + "\t" + sample + "\t" + str(c) + "\t" + allele + "\n")
allele_out.close()