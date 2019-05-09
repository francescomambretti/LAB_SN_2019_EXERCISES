CC = g++
CFLAGS = -Wall -O3 --std=c++11

ising1D.x : ising1D.o functions.o
	$(CC) functions.o ising1D.o -o ising1D.x

ising1D.o : ising1D.cpp
	$(CC) -c ising1D.cpp -o ising1D.o $(CFLAGS)

functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

clean :
	rm *.o ising1D.x
