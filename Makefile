CC = gcc
CFLAGS = -O3 -march=native -ffast-math
DEPS = tsc_x86.h
TARGET = generalized_floyd_warshall.c
OBJ = generalized_floyd_warshall.o 

%.o: $(TARGET) $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ffw: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o

bench:
	python3 benchmark.py

min_plus: clean ffw
	./ffw 256 0