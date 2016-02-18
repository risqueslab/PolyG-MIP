import sys
import os.path

notes_counter = 0
output_file = open("notes.combined.txt", 'w')

#output header line with filenames
output_file.write ("Sample\tFiltered\tSoftclipping\tUnmapped\tOff Target\tRejected SAM Flags\n")

for filename in sys.argv:
    if notes_counter > 0:
        output_file.write (os.path.basename(filename).split("_")[0] + "\t",)
        notes_file = open(filename)
        for line in notes_file:
            output_file.write (line.split(" ")[0] + "\t",)
        output_file.write("\n")
        notes_file.close()
    notes_counter += 1

output_file.close()