CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: main_8_1.exe main_8_2.exe main_8_3.exe

main_8_1.exe: main_8_1.o random.o
	$(CC) random.o main_8_1.o -o main_8_1.exe

main_8_2.exe: main_8_2.o random.o
	$(CC) random.o main_8_2.o -o main_8_2.exe

main_8_3.exe: main_8_3.o random.o
	$(CC) random.o main_8_3.o -o main_8_3.exe

main_8_1.o: main_8_1.cpp
	$(CC) -c main_8_1.cpp -o main_8_1.o $(CFLAGS)

main_8_2.o: main_8_2.cpp
	$(CC) -c main_8_2.cpp -o main_8_2.o $(CFLAGS)

main_8_3.o: main_8_3.cpp
	$(CC) -c main_8_3.cpp -o main_8_3.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -f *.o *.exe seed.out

