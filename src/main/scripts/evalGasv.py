#!/usr/bin/env python

import sys
import subprocess

gasv_filename = sys.argv[1]
truth_filename = sys.argv[2]

score_values = []

gasv_file = open(gasv_filename, "r")
for line in gasv_file:
    if line.startswith("#"):
        continue
    fields = line.split("\t")
    call_type = fields[7]
    if call_type != "D":
        continue

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

print "\t".join(["Thresh", "Calls", "TP", "Long", "WrongType", "TPR"])
for v in unique_score_values:
    calls_gte_threshold = []
    gasv_file = open(gasv_filename, "r")
    long_calls = 0
    non_del_calls = 0
    for line in gasv_file:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        call_type = fields[7]
        if float(fields[5]) >= v:
            if call_type != "D":
                non_del_calls += 1
                continue
            chrom = fields[1]
            ostart = fields[2].split(",")[0]
            oend = fields[4].split(",")[1]
            sv_len = int(oend) - int(ostart)
            if sv_len > 25000:
                long_calls += 1
                continue
            bed_line = "\t".join([chrom, ostart, oend]) + "\n"
            #print bed_line.strip()
            calls_gte_threshold.append(bed_line)
    bedtoolsProcess = subprocess.Popen(["intersectBed", "-b", "stdin", "-a", truth_filename, "-u", "-wa", "-f", ".4"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    bed_lines = ""
    for hline in calls_gte_threshold:        
        bed_lines += hline
    pout = bedtoolsProcess.communicate(bed_lines)[0]
    matches = 0
    for line in pout.split("\n"):        
        #print line
        matches += 1
    bedtoolsProcess.stdout.close()
    tpr = float(matches) / (long_calls + len(calls_gte_threshold))
    print "\t".join(map(str, [v, long_calls + len(calls_gte_threshold), matches, long_calls, non_del_calls, tpr]))
    
    
