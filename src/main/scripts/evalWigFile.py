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

def open_file(wig_filename):
    if (wig_filename.endswith("gz")):
        wig_file = gzip.open(wig_filename, "rb")
    else:
        wig_file = open(wig_filename, "r")
    return wig_file

# got the idea for this from http://www.biostars.org/post/show/6544/selecting-random-pairs-from-fastq/

records = sum(1 for _ in open_file(wig_filename))
#records = 1000
sample_size = 800000

sys.stderr.write("Total of " + str(records) + " records.\n")

rand_records = sorted([random.randint(0, records - 1) for _ in xrange(sample_size)])

wig_file = open_file(wig_filename)
rec_num = -1

sample = [0] * sample_size
sample_num = 0
for rr in rand_records:
    while rec_num < rr:
        rec_num += 1
        wig_file.readline()
    line = wig_file.readline()
    while line.startswith("variable"):
        line = wig_file.readline()
    val = float(line.split()[1])
    #print "Assigning " + str(sample_num) + " to " + str(val)
    sample[sample_num] = val
    rec_num += 1
    sample_num += 1

sample.sort()

num_quantiles = 400
quantiles = [0] * num_quantiles
q_num = 0

for i in xrange(len(sample)):
    if i % (sample_size / num_quantiles) == 0:
        quantiles[q_num] = sample[i]
        q_num += 1
    
sys.stderr.write(str(quantiles))
sys.stderr.write("\n")


def process_quantile(q):
    sys.stderr.write("quantile " + str(q) + "\n")
    temp_file_name = "/l2/users/whelanch/gene_rearrange/svpipeline/NA07051/tmp/tmp_" + str(q) + ".bed"
    extract_regions_cmd = ['condor_run', 'hadoop', 'jar', '/l2/users/whelanch/gene_rearrange/svpipeline/lib/svpipeline-1.0-SNAPSHOT-exe.jar', 'extractPositiveRegionsFromWig', '--inputWigFile', wig_filename, '--outputBedFile', temp_file_name, '--name', "tmp_" + str(q), "--faidx", faidx_filename, "--threshold", str(q), "--medianFilterWindow", medianFilterWindow]    
    subprocess.call(extract_regions_cmd)
    
    num_predictions = 0
    predicted_region = 0

    for line in open_file(temp_file_name):
        if line.startswith("track"):
            continue        
        fields = line.split()
        num_predictions += 1
        predicted_region += int(fields[2]) - int(fields[1])
    
    sys.stderr.write("num_predictions = " + str(num_predictions) + "\n")
    sys.stderr.write("predicted_region = " + str(predicted_region) + "\n")
    
    compare_bed_to_truth_cmd = ['condor_run', 'intersectBed', '-b', temp_file_name, '-a', truth_filename, '-wa', '-u']
    
    num_matches = sum(1 for _ in subprocess.Popen(compare_bed_to_truth_cmd, stdout=subprocess.PIPE).stdout)
    
    sys.stderr.write("matches = " + str(num_matches) + "\n")
    os.remove(temp_file_name)
    return (q, num_predictions, predicted_region, num_matches)

p=Pool(50)
results = p.map(process_quantile, quantiles)

for q in results:
    print "\t".join(map(str, q))
