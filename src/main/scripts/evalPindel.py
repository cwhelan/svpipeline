#!/usr/bin/env python

import sys
import subprocess
import evalBedFile

# expects bed input - chr, start, stop, and the S1 score from Pindel

pindel_filename = sys.argv[1]
truth_filename = sys.argv[2]

score_values = []

pindel_file = open(pindel_filename, "r")
for line in pindel_file:
    fields = line.split("\t")

    score = float(fields[3])
    score_values.append(score)

pindel_file.close()

unique_score_values = list(set(score_values))
unique_score_values.sort()

print "\t".join(["Thresh", "Calls", "TP", "Long", "WrongType", "TPR"])
for v in unique_score_values:
    calls_gte_threshold = []
    pindel_file = open(pindel_filename, "r")
    long_calls = 0
    non_del_calls = 0
    for line in pindel_file:
        fields = line.split("\t")
        if float(fields[3]) >= v:
            chrom = fields[0]
            ostart = fields[1]
            oend = fields[2]
            sv_len = int(oend) - int(ostart)
            if sv_len > 25000:
                long_calls += 1
                continue
            bed_line = "\t".join([chrom, ostart, oend])
            #print bed_line.strip()
            calls_gte_threshold.append(bed_line)
    (qualified_calls, matches, short_calls) = evalBedFile.eval_bed(truth_filename, calls_gte_threshold)
    tpr = float(matches) / (long_calls + len(qualified_calls))
    print "\t".join(map(str, [v, long_calls + len(qualified_calls), matches, long_calls, non_del_calls, short_calls, tpr]))
    
    
