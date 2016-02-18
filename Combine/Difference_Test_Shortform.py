#Written By: Alexander M. West
#Version: Sept. 8, 2015

#Useage: this ShortformCombinedFile Test_List #OfExperiments

import sys

output_file = open("shortform.differencetested.txt", 'w')

shortform_file = open(sys.argv[1])
current_mip = ""
samples = {}
counts = {}
minimum_reads = 40
experiments = int(sys.argv[3])
d25 = [0] * experiments
d20 = [0] * experiments
d15 = [0] * experiments
total = [0] * experiments

def output_mip():
    #print(mip)
    #print(counts)
    #print(samples)
    test_list = open(sys.argv[2])
    for line in test_list:
        number, sample1, sample2 = line.rstrip().split()
        greatest_difference = 0
        count_total1 = 0
        count_total2 = 0
        found = False
        for i in range (4,24):
            if (sample1, i) in samples and (sample2, i) in samples:
                found = True
                d = abs(samples[sample1, i] - samples[sample2, i])
                if d > greatest_difference:
                    greatest_difference = d
            elif (sample1, i) in samples:
                found = True
                d = samples[sample1, i]
                if d > greatest_difference:
                    greatest_difference = d
            elif (sample2, i) in samples:
                found = True
                d = samples[sample2, i]
                if d > greatest_difference:
                    greatest_difference = d
            if (sample1, i) in counts:
                count_total1 += counts[(sample1, i)]
            if (sample2, i) in counts:
                count_total2 += counts[(sample2, i)]
        if not found:
            output_file.write(sample1 + " vs. " + sample2 + " - sample not found\n")
        elif count_total1 < minimum_reads or count_total2 < minimum_reads:
            output_file.write(sample1 + " vs. " + sample2 + " - insufficient tags\n")
        else:
            output_file.write(sample1 + " vs. " + sample2 + " - greatest difference: " + str(greatest_difference))
            total[int(number)] += 1
            if greatest_difference < .15:
                output_file.write(" NULL\n")
            else:
                output_file.write(" MUTANT?\n")
                if greatest_difference >= .25: d25[int(number)] += 1
                if greatest_difference >= .20: d20[int(number)] += 1
                if greatest_difference >= .15: d15[int(number)] += 1
    test_list.close()

for line in shortform_file:
    mip, sample, allele, length, tags, proportion = line.rstrip().split()
    if mip != current_mip:
        output_mip()
        samples = {}
        counts = {}
        current_mip = mip
        output_file.write(current_mip +"\n")
    if (sample, int(length)) not in samples:
        samples[(sample, int(length))] = float(proportion)
        counts[(sample, int(length))] = int(tags)
    else:
        samples[(sample, int(length))] += float(proportion)
        counts[(sample, int(length))] += int(tags)
output_mip()

output_file.write("\n\n\n*** Summary ***\n")
output_file.write("d > .15 \t" + str(d15) + "\n")
output_file.write("d > .20 \t" + str(d20) + "\n")
output_file.write("d > .25 \t" + str(d25) + "\n")
output_file.write("total \t" + str(total) + "\n")

shortform_file.close()
output_file.close()
