CC        = cc
CFLAGS    = -Wall -O3 -march=native
NAME      = ucrdtw

all: bin/ucrdtw

bin/ucrdtw: ucrdtw.o
	$(CC) $(CFLAGS) ucrdtw.o -o bin/ucrdtw

ucrdtw.o: src/ucrdtw.c src/ucrdtw.h
	$(CC) $(CFLAGS) -c src/ucrdtw.c 

clean:
	rm -f ucrdtw.o bin/ucrdtw