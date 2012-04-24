#!/usr/bin/env python

import sys
import subprocess

breakdancer_filename = sys.argv[1]
truth_filename = sys.argv[2]

score_values = []

breakdancer_file = open(breakdancer_filename, "r")
for line in breakdancer_file:
    fields = line.split("\t")
    final_weighted_score = float(fields[7])
    score_values.append(final_weighted_score)
breakdancer_file.close()

unique_score_values = list(set(score_values))
unique_score_values.sort()

print "\t".join(["Thresh", "Calls", "TP", "Long", "TPR"])
for v in unique_score_values:
    calls_gte_threshold = []
    breakdancer_file = open(breakdancer_filename, "r")
    long_calls = 0
    for line in breakdancer_file:
        if float(line.split("\t")[7]) >= v:
            sv_len = int(line.split("\t")[5]) - int(line.split("\t")[1])
            if sv_len > 10000:
                long_calls += 1
                continue
            calls_gte_threshold.append(line)
    bedtoolsProcess = subprocess.Popen(["pairToBed", "-type", "ispan",  "-a", "stdin", "-b", truth_filename], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for hline in calls_gte_threshold:
        bedtoolsProcess.stdin.write(hline)
    bedtoolsProcess.stdin.close()
    matches = 0
    for line in bedtoolsProcess.stdout:        
        #print line
        matches += 1
    bedtoolsProcess.stdout.close()
    print "\t".join(map(str, [v, long_calls + len(calls_gte_threshold), matches, long_calls, float(matches) / (long_calls + len(calls_gte_threshold))]))
    
    
