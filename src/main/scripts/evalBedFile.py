#!/usr/bin/env python

import sys
import subprocess

def eval_bed(truth_filename, predictions):
    bedtools_process = subprocess.Popen(["intersectBed", "-a", "stdin", "-b", truth_filename, "-u"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    pstdout = bedtools_process.communicate("\n".join(predictions))[0]
    matches = 0
    for line in pstdout.split("\n"):
        #	print "line: " + line
        matches += 1
        # for the last newline
    matches = matches - 1
    return matches

# this script evaluates a series of bed lines from stdin against a truth file (passed as an argument) and returns the number of correct calls
if __name__ == "__main__":
    import sys

    # slurp all the bed lines into memory
    predictions = []
    for line in sys.stdin:
        predictions.append(line)

    eval_bed(sys.argv[1], predictions)
