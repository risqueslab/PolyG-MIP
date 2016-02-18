import sys
import os.path

significant_figures = 5
shortform_counter = 0
combined_allele_count = {} #mip,allele tuple key; proportion list value
output_file = open("shortform.combined.txt", 'w')

#output header line with filenames
output_file.write ("MIP\tallele\tlengtif h\t",)
for filename in sys.argv:
    if shortform_counter > 0:
        output_file.write (os.path.basename(filename).partition("_")[0] + " Tags\t",)
        output_file.write (os.path.basename(filename).partition("_")[0] + " %\t",)
    shortform_counter += 1
output_file.write("\n")
shortform_counter = 0

for filename in sys.argv:
    if shortform_counter > 0:
        shortform_file = open(filename)
        for line in shortform_file:
            mip, barcode, allele, length, reads, proportion = line.rstrip().split()
            if (mip, allele, length) in combined_allele_count:
                combined_allele_count[mip, allele, length][shortform_counter - 1] = (reads, proportion)
            else:
                combined_allele_count[mip, allele, length] = [(0,0)] * (len(sys.argv) -1)
                combined_allele_count[mip, allele, length][shortform_counter - 1] = (reads, proportion)
        shortform_file.close()
    shortform_counter += 1

for mip_allele in combined_allele_count:
    output_file.write (str(mip_allele[0]) + "\t" + str(mip_allele[1]) + "\t" + str(mip_allele[2]) + "\t",)
    for (reads, proportion) in (combined_allele_count[mip_allele]):
        if (reads, proportion) == (0,0):
            output_file.write(" \t \t",)
        else:
            #output_file.write(str(item) + "\t",)
            output_file.write(str(reads) + "\t" + str(round(float(proportion), significant_figures)) + "\t",)
    output_file.write("\n")
