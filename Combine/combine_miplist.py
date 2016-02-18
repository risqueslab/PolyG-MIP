import sys
import os.path

miplist_counter = 0
combined_mip_data = {} #mip key, read/unused tag/used tag tuple data
output_file = open("miplist.combined.txt", 'w')

#output header line with filenames
output_file.write ("MIP\t",)
for filename in sys.argv:
    if miplist_counter > 0:
        #output_file.write (os.path.basename(filename).split("_")[0] + " Reads\t",)
        output_file.write (os.path.basename(filename).split(".")[2] + " Reads\t",)
        output_file.write (os.path.basename(filename).split(".")[2] + " Unused Tags\t",)
        output_file.write (os.path.basename(filename).split(".")[2] + " Used Tags\t",)
        #output_file.write (os.path.basename(filename).split("_")[0] + " Tags < 3\t",)
        #output_file.write (os.path.basename(filename).split("_")[0] + " Tags > 2\t",)
    miplist_counter += 1
output_file.write("\n")
miplist_counter = 0

for filename in sys.argv:
    if miplist_counter > 0:
        miplist_file = open(filename)
        next(miplist_file)
        for line in miplist_file:
            mip, reads, unused_tags, used_tags = line.rstrip().split()
            if mip in combined_mip_data:
                combined_mip_data[mip][miplist_counter - 1] = (reads, unused_tags, used_tags)
            else:
                combined_mip_data[mip] = [(0, 0, 0)] * (len(sys.argv) -1)
                combined_mip_data[mip][miplist_counter - 1] = (reads, unused_tags, used_tags)
        miplist_file.close()
    miplist_counter += 1

for mip in combined_mip_data:
    output_file.write (str(mip) + "\t",)
    for (reads, unused_tags, used_tags) in combined_mip_data[mip]:
        output_file.write(str(reads) + "\t" + str(unused_tags) + "\t" + str(used_tags) + "\t",)
    output_file.write("\n")