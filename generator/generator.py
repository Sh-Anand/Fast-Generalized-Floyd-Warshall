from constants import program

f = open("generated_vectorized_tiled.c", "w");

f.write(program)

f.close()