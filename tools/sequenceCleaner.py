#sequence cleaner removes sequences that are ambiguous (6-mer appending the poly sequence is indefinite ("N") and shifts all "N" characters
#in poly sequence right so that they can be combined

def sequenceCleaner(string):
	if len(string) < 13:
		return "", 0
	if string[5] == "*":
		return "", 0
	if string[len(string)-6] == "*":
		return "", 0
	if not "*" in string[6:len(string)-6]:
		return string[6:len(string)-6], len(string)-12
	else:
		return rightShiftN(string[6:len(string)-6])
		
def rightShiftN(string):
	output = ""
	indefinite = ""
	total = 0.0
	for i in string:
		if i == "*":
			indefinite = indefinite + i
			total += 0.51
		else:
			output = output + i
			total += 1.0
	output = output + indefinite
	return output, total