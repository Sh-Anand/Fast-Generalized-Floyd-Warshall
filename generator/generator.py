import sys
from subprocess import run, PIPE
import writer
import writer_or_and
import matplotlib.pyplot as plt

config_file = "../best_result_dump_new.csv"
generated_file = "generated_vectorized_tiled.c"

op = int(sys.argv[1])
if len(sys.argv) == 3:
    f = open(generated_file, "w")
    program = ""
    if op == 2:
        program = writer_or_and.generate_program(int(sys.argv[2]), 0)
    else:
        program = writer.generate_program(int(sys.argv[2]), op)
    f.write(program)
    f.close()
elif len(sys.argv) == 2:
    max_speedup = 0
    best_unrolling = 0
    repetitions_for_confidence = 3
    x = []
    y = []
    for unroll in (2 ** p for p in range(0, 5)):
        x.append(unroll)
        f = open(generated_file, "w")
        program = ""
        if op == 2:
            program = writer_or_and.generate_program(unroll, 0)
        else:
            program = writer.generate_program(unroll, op)
        f.write(program)
        f.close()
        run("gcc -O3 -ffast-math -march=native " + generated_file, shell=True)

        num = 0
        mean_speedup = 0
        fwi = -1
        with open(config_file, "r") as config:
            lines = config.readlines()
            for line in lines:
                # Handle each semi-ring
                if len(line.split(",")) < 2:
                    fwi = fwi + 1

                # Find corresponding semi-ring data in config file
                if op == fwi:
                    (n, l1, b, _) = line.split(",")

                    speedup = 0
                    for reps in range(repetitions_for_confidence):
                        out = run("./a.out "+n+" "+l1+" "+b, shell=True, stdout=PIPE).stdout.decode('utf-8')
                        vals = [float(t) for t in out.split()]
                        baseline_t = vals[0]
                        generated_t = vals[1]
                        speedup = speedup + baseline_t/generated_t
                    mean_speedup = mean_speedup + (speedup/repetitions_for_confidence)
                    num = num + 1

        mean_speedup = mean_speedup/num
        y.append(mean_speedup)
        if mean_speedup > max_speedup:
            max_speedup = mean_speedup
            best_unrolling = unroll

    plt.plot(x, y)
    plt.xlabel("Unrolling")
    plt.ylabel("Speedup")
    plt.savefig("unrolling_benchmark.png")
    plt.show()
    print("Best unrolling = " + str(best_unrolling) + " with speedup = " + str(max_speedup))
else:
    print("Supply exactly one argument: unrolling factor")
