#!/usr/bin/env python

import sys
import subprocess

gasv_filename = sys.argv[1]
truth_filename = sys.argv[2]

score_values = []

gasv_file = open(gasv_filename, "r")
for line in gasv_file:
    fields = line.split("\t")
    final_weighted_score = float(fields[5])
    score_values.append(final_weighted_score)
    # chrom = fields[1]
    # ostart = fields[2].split(",")[0]
    # oend = fields[4].split(",")[1]            
    # bed_line = "\t".join([chrom, ostart, oend])
    # print bed_line

gasv_file.close()

unique_score_values = list(set(score_values))
unique_score_values.sort()

print "\t".join(["Thresh", "Calls", "TP", "Long", "TPR"])
for v in unique_score_values:
    calls_gte_threshold = []
    gasv_file = open(gasv_filename, "r")
    long_calls = 0
    for line in gasv_file:
        fields = line.split("\t")
        if float(fields[5]) >= v:
            chrom = fields[1]
            ostart = fields[2].split(",")[0]
            oend = fields[4].split(",")[1]
            sv_len = int(oend) - int(ostart)
            if sv_len > 100000:
                long_calls += 1
                continue
            bed_line = "\t".join([chrom, ostart, oend]) + "\n"
            calls_gte_threshold.append(bed_line)
    bedtoolsProcess = subprocess.Popen(["intersectBed", "-b", "stdin", "-a", truth_filename, "-u", "-wa"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for hline in calls_gte_threshold:        
        bedtoolsProcess.stdin.write(hline)
    bedtoolsProcess.stdin.close()
    matches = 0
    for line in bedtoolsProcess.stdout:        
        #print line
        matches += 1
    bedtoolsProcess.stdout.close()
    print "\t".join(map(str, [v, long_calls + len(calls_gte_threshold), matches, long_calls, float(matches) / (long_calls + len(calls_gte_threshold))]))
    
    
