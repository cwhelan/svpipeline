#!/usr/bin/env python

import sys
import subprocess
import evalBedFile

hydra_filename = sys.argv[1]
truth_filename = sys.argv[2]

support_values = []

hydra_file = open(hydra_filename, "r")
for line in hydra_file:
    fields = line.split("\t")
    final_weighted_support = float(fields[18])
    support_values.append(final_weighted_support)
hydra_file.close()

unique_support_values = list(set(support_values))
unique_support_values.sort()
#print unique_support_values

print "\t".join(["Thresh", "Calls", "TP", "Long", "TPR", "WrongType"])
for v in unique_support_values:
    calls_gte_threshold = []
    hydra_file = open(hydra_filename, "r")
    long_calls = 0
    wrong_type = 0
    for line in hydra_file:
        fields = line.split("\t")
        if float(fields[18]) >= v:
#            print line
            s1 = fields[8]
            s2 = fields[9]
#            print "strands " + s1 + s2
            if not ((s1 == "+" and s2 == "-") or (s1 == "-" and s2 == "+")):
#                print "wrong type: " + line
                wrong_type += 1
                continue 
            sv_len = int(fields[14])
            if sv_len > 10000:
                long_calls += 1
                continue
            calls_gte_threshold.append(line)
    bed_lines = []
    for line in calls_gte_threshold:
        fields = line.split("\t")
        bed_line = "\t".join([fields[0], fields[2], fields[4]])
        bed_lines.append(bed_line)

    matches = evalBedFile.eval_bed(truth_filename, bed_lines)
    print "\t".join(map(str, [v, long_calls + len(calls_gte_threshold), matches, long_calls, float(matches) / (long_calls + len(calls_gte_threshold)), wrong_type]))
    
    
