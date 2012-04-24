#!/usr/bin/env python

import sys
import subprocess

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

print "\t".join(["Thresh", "Calls", "TP", "Long", "TPR"])
for v in unique_support_values:
    calls_gte_threshold = []
    hydra_file = open(hydra_filename, "r")
    long_calls = 0
    for line in hydra_file:
        if float(line.split("\t")[18]) >= v:
            sv_len = int(line.split("\t")[14])
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
    
    
