#Written By: Alexander M. West
#Version: Sept. 8, 2015

#Useage: path ShortformFile Test_List #OfExperiments

import sys
from scipy.stats import spearmanr

output_file = open("shortform.spearmantested.txt", 'w')

shortform_file = open(sys.argv[1])
current_mip = ""
samples = {}
counts = {}
experiments = int(sys.argv[3])
c05 = [0] * experiments
c03 = [0] * experiments
c01 = [0] * experiments
total = [0] * experiments

def output_mip():
    #print(mip)
    #print(counts)
    #print(samples)
    test_list = open(sys.argv[2])
    for line in test_list:
        number, sample1, sample2 = line.rstrip().split()
        array1 = []
        array2 = []
        found = 0
        for i in range (4,24):
            if (sample1, i) in samples and (sample2, i) in samples:
                found += 1
                array1 += [samples[sample1, i]]
                array2 += [samples[sample2, i]]
            elif (sample1, i) in samples:
                found += 1
                array1 += [samples[sample1, i]]
                array2 += [0]
            elif (sample2, i) in samples:
                found += 1
                array1 += [0]
                array2 += [samples[sample2, i]]
        #print(str(array1))
        #print(str(array2))
        if not found > 1:
            output_file.write(sample1 + " vs. " + sample2 + " - sample not found\n")
        #elif count_total1 < minimum_reads or count_total2 < minimum_reads:
        #output_file.write(sample1 + " vs. " + sample2 + " - insufficient tags\n")
        else:
            (c, p) = spearmanr(array1, array2)
            output_file.write(sample1 + " vs. " + sample2 + " - p-value: " + str(round(float(1.0) - p, 4)) + " - c-value: " + str(round(float(1.0) - c, 4)))
            total[int(number)] += 1
            if abs(c) < .5:
                output_file.write(" MUTANT?\n")
                if abs(c) < .5: c05[int(number)] += 1
                if abs(c) < .3: c03[int(number)] += 1
                if abs(c) < .1: c01[int(number)] += 1
            else:
                output_file.write("\n")
    test_list.close()

for line in shortform_file:
    mip, sample, allele, length, tags, proportion = line.rstrip().split()
    if mip != current_mip:
        output_mip()
        samples = {}
        current_mip = mip
        output_file.write(current_mip +"\n")
    if (sample, int(length)) not in samples:
        samples[(sample, int(length))] = int(tags)
    else:
        samples[(sample, int(length))] += int(tags)
output_mip()

output_file.write("\n\n\n*** Summary ***\n")
output_file.write("c < .5 \t" + str(c05) + "\n")
output_file.write("c < .3 \t" + str(c03) + "\n")
output_file.write("c < .1 \t" + str(c01) + "\n")
output_file.write("totals \t" + str(total) + "\n")

shortform_file.close()
output_file.close()
