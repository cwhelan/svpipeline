#!/usr/bin/env python

import sys
import subprocess

breakdancer_filename = sys.argv[1]
truth_filename = sys.argv[2]
#bed_file = sys.argv[3]

score_values = []

breakdancer_file = open(breakdancer_filename, "r")
for line in breakdancer_file:
    if line.startswith("#"):
        continue
    fields = line.split("\t")
    score = float(fields[8])
    score_values.append(score)
breakdancer_file.close()

unique_score_values = list(set(score_values))
unique_score_values.sort()

print "\t".join(["Thresh", "Calls", "TP", "Long", "WrongType", "TPR"])
for v in unique_score_values:
    calls_gte_threshold = []
    breakdancer_file = open(breakdancer_filename, "r")
    long_calls = 0
    non_del_calls = 0
    for line in breakdancer_file:
#        print line
        fields = line.split("\t")
        if line.startswith("#"):
#            print "comment!"
            continue
        if float(fields[8]) >= v:
#            print "gte v"
            if (not fields[6] == "DEL") or (fields[0] != fields[3]):
                non_del_calls += 1
                continue
            sv_len = int(fields[7])
            # sys.stderr.write("len: " + str(sv_len) + "\n")
            if sv_len > 25000:
#                print "too long"
                long_calls += 1
                continue
            calls_gte_threshold.append(line)
#    sys.stderr.write(str(calls_gte_threshold))
    bedtoolsProcess = subprocess.Popen(["intersectBed", "-a", "stdin", "-b", truth_filename], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    for line in calls_gte_threshold:
        fields = line.split("\t")
        bed_line = "\t".join([fields[0], fields[1], fields[4]]) + "\n"        
        bedtoolsProcess.stdin.write(line)
    bedtoolsProcess.stdin.close()
    matches = 0
    for line in bedtoolsProcess.stdout:        
        #print line
        matches += 1
    bedtoolsProcess.stdout.close()
    # sys.stderr.write("v = " + str(v))
    # sys.stderr.write("long_calls = " + str(long_calls))
    # sys.stderr.write("calls_gte_thresh = " + str(calls_gte_threshold))
    calls = long_calls + len(calls_gte_threshold)
    tpr = float(matches) / calls
    print "\t".join(map(str, [v, calls, matches, long_calls, non_del_calls, tpr]))
    
    
