import sys
from writer import generate_program

if len(sys.argv) == 2:
    f = open("generated_vectorized_tiled.c", "w");
    program = generate_program(int(sys.argv[1]))
    f.write(program)
    f.close()
else:
    print("Supply exactly one argument: unrolling factor")