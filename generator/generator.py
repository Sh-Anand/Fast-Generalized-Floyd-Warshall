import string
import sys
from subprocess import run, PIPE
import writer
import writer_or_and
import matplotlib.pyplot as plt

config_file = "best_result_dump_new.csv"
fw = ["min_plus", "max_min", "or_and"]


def write_generated_code(file_name: string, op_idx: int, unroll_size: int):
    f = open(file_name, "w")
    program = ""
    if op_idx == 2:
        program = writer_or_and.generate_program(unroll_size, 0)
    else:
        program = writer.generate_program(unroll_size, op_idx)
    f.write(program)
    f.close()


op = int(sys.argv[1])
generated_file = fw[op]+"_generated_vectorized_tiled.c"

if len(sys.argv) == 3:
    write_generated_code(generated_file, op, int(sys.argv[2]))
elif len(sys.argv) == 2:
    max_speedup = 0
    best_unrolling = 0
    repetitions_for_confidence = 3
    x = []
    y = []
    for unroll in range(1, 13):
        try:
            x.append(unroll)
            write_generated_code(generated_file, op, unroll)
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
                        continue

                    # Find corresponding semi-ring data in config file
                    if op == fwi:
                        (n, l1, b, _) = line.split(",")
                        if int(n) < 216:
                            continue

                        speedup = 0
                        for reps in range(repetitions_for_confidence):
                            out = run("./a.out "+n+" "+l1+" "+b, shell=True, stdout=PIPE).stdout.decode('utf-8')
                            vals = [float(t) for t in out.split()]
                            baseline_t = vals[1]
                            generated_t = vals[0]
                            speedup = speedup + baseline_t/generated_t
                        mean_speedup = mean_speedup + (speedup/repetitions_for_confidence)
                        num = num + 1

            mean_speedup = mean_speedup/num
            y.append(mean_speedup)
            if mean_speedup > max_speedup:
                max_speedup = mean_speedup
                best_unrolling = unroll

        except IndexError:
            continue

    # Generate code with best found unrolling factor
    write_generated_code(generated_file, op, best_unrolling)

    plt.plot(x, y)
    plt.xlabel("Unrolling")
    plt.ylabel("Speedup")
    plt.savefig("unrolling_benchmark.png")
    plt.show()
    print("Best unrolling for " + fw[op] + " is = " + str(best_unrolling) + " with speedup = " + str(max_speedup))
else:
    print("Supply exactly one argument: unrolling factor")
