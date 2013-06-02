#!/usr/bin/env python

import sys
import subprocess
import evalBedFile

# Arguments:
#
# Breakdancer file: the bd_out.txt file generated by breakdancer
# Truth file: File with the true deletions

breakdancer_filename = sys.argv[1]
truth_filename = sys.argv[2]
#bed_file = sys.argv[3]

score_values = []

breakdancer_file = open(breakdancer_filename, "r")
for line in breakdancer_file:
    if line.startswith("#"):
        continue
    fields = line.split("\t")
    score = float(fields[9])
    score_values.append(score)
breakdancer_file.close()

unique_score_values = list(set(score_values))
unique_score_values.sort()

print "\t".join(["Thresh", "Calls", "TP", "WrongType", "Short", "TPR"])
for v in unique_score_values:
    calls_gte_threshold = []
    breakdancer_file = open(breakdancer_filename, "r")
    non_del_calls = 0
    # sometimes Breakdancer makes nonsensical calls; deletions with ends on different chromosomes, we have to
    # track those separately here
    bad_calls = 0
    for line in breakdancer_file:
#        print line
        fields = line.split("\t")
        if line.startswith("#"):
#            print "comment!"
            continue
        if float(fields[9]) >= v:
#            print "gte v"
            if not (fields[6] == "DEL"):
                non_del_calls += 1
                continue
            if (fields[0] != fields[3]):
                bad_calls += 1
                continue
            sv_len = int(fields[7])
            # sys.stderr.write("len: " + str(sv_len) + "\n")
            calls_gte_threshold.append(line)
            #    sys.stderr.write(str(calls_gte_threshold))

    bed_lines = []
    for line in calls_gte_threshold:
        fields = line.split("\t")
        bed_line = "\t".join([fields[0], fields[1], fields[4]])
        bed_lines.append(bed_line)

    (qualified_calls, matches, short_calls) = evalBedFile.eval_bed(truth_filename, bed_lines)
    qualified_calls += bad_calls
    tpr = float(matches) / (qualified_calls)
    print "\t".join(map(str, [v, qualified_calls, matches, non_del_calls, short_calls, tpr]))
    
    
