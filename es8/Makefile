CC = g++
CFLAGS = -Wall -O3 --std=c++11

Main2.exe : Main2.o random.o
	$(CC) random.o Main2.o -o Main2.exe
Main2.o : Main2.cpp
	$(CC) -c Main2.cpp -o Main2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Main2.exe
