CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : analize.o
	$(CC) analize.o -o main.exe
analize.o : analize.cpp
	$(CC) -c analize.cpp -o analize.o $(CFLAGS)
clean :
	rm *.o main.exe