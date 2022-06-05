import sys
from subprocess import run, PIPE
from csv import DictReader
from writer import generate_program
import matplotlib.pyplot as plt

if len(sys.argv) == 2:
    f = open("generated_vectorized_tiled.c", "w")
    program = generate_program(int(sys.argv[1]))
    f.write(program)
    f.close()
elif len(sys.argv) == 1:
    max_speedup = 0
    best_unrolling = 0
    maxreps = 3
    for unroll in (2 ** p for p in range(0, 5)):
        f = open("generated_vectorized_tiled.c", "w")
        program = generate_program(unroll)
        f.write(program)
        f.close()
        run("gcc -O3 -ffast-math -march=native generated_vectorized_tiled.c", shell=True)

        test_file = open("generator_args.csv", "r")
        tests = DictReader(test_file)
        num = 0
        mean_speedup = 0
        for test in tests:
            speedup = 0
            for reps in range(maxreps):
                out = run("./a.out "+str(test['N'])+" "+str(test['L1'])+" "+str(test['B']), shell=True, stdout=PIPE).stdout.decode('utf-8')
                vals = [float(t) for t in out.split()]
                print(vals)
                speedup = speedup + vals[0]/vals[1]
            mean_speedup = mean_speedup + speedup/maxreps
            num = num + 1
        
        mean_speedup = mean_speedup/num
        if mean_speedup > max_speedup:
            mean_speedup = max_speedup
            best_unrolling = unroll
        test_file.close()
        print(str(unroll) + " Done")
    
    print("Best unrolling = " + str(best_unrolling) + " with speedup = " + str(speedup))
else:
    print("Supply exactly one argument: unrolling factor")