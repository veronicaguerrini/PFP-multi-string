# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c11 -g
CC=gcc

# executables not using threads (and therefore not needing the thread library)
EXECS_NT=newscan.x bwtparse bwtparse64 pfmultistringbwt.x pfmultistringbwt64.x

# targets not producing a file declared phony
.PHONY: all clean lcp

all: $(EXECS_NT)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

bwtparse: bwtparse.c gsa/gsacak.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

bwtparse64: bwtparse.c gsa/gsacak64.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl -DM64

newscan.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

# prefix free multi-string BWT construction
pfmultistringbwt.x: pfmultistringbwt.cpp gsa/gsacak.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

# as above for large files
pfmultistringbwt64.x: pfmultistringbwt.cpp gsa/gsacak64.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DM64

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS_NT) *.o gsa/*.o

