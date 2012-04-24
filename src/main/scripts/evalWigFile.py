#!/usr/bin/env python

import sys
import gzip
import random
import subprocess
import tempfile
import os
from multiprocessing import Pool
import sys

wig_filename = sys.argv[1]
truth_filename = sys.argv[2]
faidx_filename = sys.argv[3]
medianFilterWindow = sys.argv[4]
lower_threshold = float(sys.argv[5])

def open_file(wig_filename):
    if (wig_filename.endswith("gz")):
        wig_file = gzip.open(wig_filename, "rb")
    else:
        wig_file = open(wig_filename, "r")
    return wig_file

values_above_threshold = []
wig_file = open_file(wig_filename)
for line in wig_file:
    if line.startswith("track") or line.startswith("variable"):
        continue
    val = float(line.split()[1])
    if val > lower_threshold:
        values_above_threshold.append(val)

#print "values above threshold: " + str(len(values_above_threshold))
values_above_threshold.sort()

num_quantiles = 500
quantiles = [0] * (num_quantiles + 1)
q_num = 0

for i in xrange(len(values_above_threshold)):    
    if i % (len(values_above_threshold) / num_quantiles) == 0:
        if (q_num > num_quantiles):
            continue
#        sys.stderr.write("i: " + str(i))
        qv = values_above_threshold[i]
#        sys.stderr.write("q_num: " + str(q_num) + " = " + str(qv))
        quantiles[q_num] = qv
        #print "quantiles " + str(q_num) + " = " + str(values_above_threshold[i])
        q_num += 1
    
sys.stderr.write(str(quantiles))
sys.stderr.write("\n")

def process_quantile(q):
    eval_at_q_cmd = ['condor_run', 'python', '/l2/users/whelanch/gene_rearrange/svpipeline/src/main/scripts/evalWigFileAtThreshold.py', str(q), wig_filename, truth_filename, faidx_filename, medianFilterWindow]
    #print eval_at_q_cmd
    result = subprocess.Popen(eval_at_q_cmd, stdout=subprocess.PIPE).communicate()[0]
    result_fields = result.split()
    num_predictions = int(result_fields[1])
    predicted_region = int(result_fields[2])
    num_matches = int(result_fields[3])
    return (q, num_predictions, predicted_region, num_matches, float(num_matches) / num_predictions)
    
p=Pool(50)
results = p.map(process_quantile, quantiles)

print "\t".join(["Thresh", "Calls", "Region", "TP", "TPR"])
for q in results:
    print "\t".join(map(str, q))
