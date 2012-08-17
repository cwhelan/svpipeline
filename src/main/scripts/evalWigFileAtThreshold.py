#!/usr/bin/env python

import sys
from cStringIO import StringIO
import random
import subprocess
import tempfile
import os
import sys

def open_file(wig_filename):
    if (wig_filename.endswith("gz")):
        p = subprocess.Popen(["zcat",wig_filename], 
                             stdout = subprocess.PIPE)
        wig_file = p.stdout
    else:
        wig_file = open(wig_filename, "r")
    return wig_file

q = float(sys.argv[1])
wig_filename = sys.argv[2]
truth_filename = sys.argv[3]
faidx_filename = sys.argv[4]
median_filter_window = sys.argv[5]

#sys.stderr.write("quantile " + str(q) + "\n")
#temp_file_name = "tmp/tmp_" + str(q) + ".bed"
temp_file = tempfile.NamedTemporaryFile()
temp_file_name = temp_file.name
extract_regions_cmd = ['hadoop', 'jar', '/l2/users/whelanch/gene_rearrange/svpipeline/lib/cloudbreak-1.0-SNAPSHOT-exe.jar', 'extractPositiveRegionsFromWig', '--inputWigFile', wig_filename, '--outputBedFile', temp_file_name, '--name', "tmp_" + str(q), "--faidx", faidx_filename, "--threshold", str(q), "--medianFilterWindow", median_filter_window]
subprocess.call(extract_regions_cmd)

num_predictions = 0
predicted_region = 0

for line in open_file(temp_file_name):
    if line.startswith("track"):
        continue        
    fields = line.split()
    num_predictions += 1
    predicted_region += int(fields[2]) - int(fields[1])
    
#sys.stderr.write("num_predictions = " + str(num_predictions) + "\n")
#sys.stderr.write("predicted_region = " + str(predicted_region) + "\n")
    
compare_bed_to_truth_cmd = ['intersectBed', '-b', temp_file_name, '-a', truth_filename, '-wa', '-u', '-f', '.4']
    
num_matches = sum(1 for _ in subprocess.Popen(compare_bed_to_truth_cmd, stdout=subprocess.PIPE).stdout)
tpr = float(num_matches) / num_predictions
    
#sys.stderr.write("matches = " + str(num_matches) + "\n")
temp_file.close()
print "\t".join(map(str, [q, num_predictions, predicted_region, num_matches, tpr]))
