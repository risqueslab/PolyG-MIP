import sys
import os.path

significant_figures = 5
shortform_counter = 0
combined_allele_count = {} #mip,allele tuple key; proportion list value
output_file = open("shortform.combined.txt", 'w')

#output header
output_file.write ("MIP\tSample\tAllele\tLength\tTags\t%\n",)

for filename in sys.argv:
    if shortform_counter > 0:
        shortform_file = open(filename)
        for line in shortform_file:
            mip, barcode, allele, length, tags, proportion = line.rstrip().split()
            output_file.write(mip + "\t" + os.path.basename(filename).partition("_")[0] + "\t" + str(allele) + "\t" + str(length) + "\t" + str(tags) + "\t" + str(round(float(proportion), significant_figures)) + "\n")
        shortform_file.close()
    shortform_counter += 1

output_file.close()
