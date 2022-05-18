from subprocess import run

import re
import os
import matplotlib.pyplot as plt
from platform import system
from statistics import median

try:
    os.makedirs("./plots")
except Exception as err:
    pass #probably exists

compiled_file = "./a" + (".exe " if system() == "Windows" else ".out ") #do not remove space

start, end, step = 10,250,10 #from n = start, to n = end, in steps of step
repetitions_for_confidence = 10 #number of repetitions of each n, median is taken
file_name = "generalized_floyd_warshall.c " #do not remove space after the end!

print(os.getcwd())
flags = ["-O0", "-O3 -fno-tree-vectorize", "-O3 -ffast-math -march=native"]
fw = ["min_plus", "or_and", "max_min"] #functions in generalized floyd warshall

for fwi in range(len(fw)):
    print("Benchmarking " + str(fw[fwi]))
    plt.clf()
    for flag in flags:
        run("gcc -mavx -w " + file_name + flag, shell=True)
        x = []
        y = []
        for n in range(start, end, step):
            flops_per_cycle = []
            for i in range(repetitions_for_confidence):
        
                output = run("./a.exe " +str(n)+" "+str(fwi), capture_output=True).stdout.decode("utf-8")

                cycles = [float(t) for t in output.split() if re.match(r'^-?\d+(?:\.\d+)$', t) is not None][0]
                # print(cycles)
                flops_per_cycle.append(2*(n**3) / cycles)

            x.append(n)
            y.append(median(flops_per_cycle))
        plt.plot(x,y,label=flag)

    print("Done benchmarking " + str(fw[fwi]))
    plt.legend()
    plt.xlabel("N")
    plt.ylabel("flops/cycle")
    plt.savefig("./plots/"+fw[fwi])
    print("Generate plots for " + str(fw[fwi]))