from subprocess import run

import re
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
from statistics import median

try:
    os.makedirs("./plots")
except Exception as err:
    pass #probably exists

start, end, step = 10,100,20 #from n = start, to n = end, in steps of step
repetitions_for_confidence = 5 #number of repetitions of each n, median is taken
file_name = "generalized_floyd_warshall.c " #do not remove space after the end!

flags = ["-O0", "-O3 -fno-tree-vectorize", "-O3 -ffast-math -march=native"]
fw = ["min_plus", "or_and", "max_min"] #functions in generalized floyd warshall

for fwi in range(len(fw)):
    plt.clf()
    for flag in flags:
        run("gcc " + file_name + flag, shell=True)
        x = []
        y = []
        for n in range(start, end, step):
            flops_per_cycle = []
            for i in range(repetitions_for_confidence):
        
                output = run("./a.exe "+str(n)+" "+str(fwi), capture_output=True).stdout.decode("utf-8")

                cycles = [float(t) for t in output.split() if re.match(r'^-?\d+(?:\.\d+)$', t) is not None][0]
                print(cycles)
                flops_per_cycle.append(2*(n**3) / cycles)

            x.append(n)
            y.append(median(flops_per_cycle))

        plt.plot(x,y,label=flag)

    plt.legend()
    plt.xlabel("N")
    plt.ylabel("flops/cycle")
    plt.savefig("./plots/"+fw[fwi])