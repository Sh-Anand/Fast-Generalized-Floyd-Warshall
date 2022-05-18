CC = gcc
CFLAGS = -O3 -march=native -ffast-math
DEPS = tsc_x86.h
TARGET = generalized_floyd_warshall.c
OBJ = generalized_floyd_warshall.o 

%.o: $(TARGET) $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ffw: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

ffw_blocked: $(OBJ)
	$(CC) -o -D name=BLOCKED $@ $^ $(CFLAGS) 

ffw_tiled: $(OBJ)
	$(CC) -D name=TILED -o $@ $^ $(CFLAGS) 

clean:
	rm -f *.o

bench:
	python3 benchmark.py

min_plus: clean ffw
	./ffw 256 0

test_blocked: clean ffw_blocked
	./ffw_blocked 256 0

test_tiled: clean ffw_tiled
	./ffw_tiled 256 0