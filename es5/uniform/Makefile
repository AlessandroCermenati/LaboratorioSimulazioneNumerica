CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : Metropolis.o random.o Autofunzione.o
	$(CC) Metropolis.o random.o Autofunzione.o -o main.exe
Metropolis.o : Metropolis.cpp random.h Autofunzione.h funzionebase.h
	$(CC) -c Metropolis.cpp -o Metropolis.o $(CFLAGS)
random.o : random.cpp random.h funzionebase.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Autofunzione.o : Autofunzione.h
	$(CC) -c Autofunzione.cpp -o Autofunzione.o $(CFLAGS)
clean :
	rm *.o main.exe *.dat
