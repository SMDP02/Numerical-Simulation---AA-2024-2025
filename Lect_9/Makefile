CC = g++
CFLAGS = -Wall -O3 --std=c++11 # o c++14 se preferisci
# Aggiungi -larmadillo se il linker non lo trova automaticamente
LIBS = -larmadillo 

# 1. Aggiunto individual.o alla lista delle dipendenze del file eseguibile
all: main_9_1.exe

main_9_1.exe: main_9_1.o random.o system.o individual.o
	$(CC) random.o system.o individual.o main_9_1.o -o main_9_1.exe $(LIBS)

main_9_1.o: main_9_1.cpp system.h random.h
	$(CC) -c main_9_1.cpp -o main_9_1.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

# 2. individual.o non ha dipendenze da system, solo dai suoi header
individual.o: individual.cpp individual.h random.h
	$(CC) -c individual.cpp -o individual.o $(CFLAGS)

# 3. system.o dipende sia dal suo header che da quello di individual
system.o: system.cpp system.h individual.h random.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)

clean:
	rm -f *.o *.exe seed.out