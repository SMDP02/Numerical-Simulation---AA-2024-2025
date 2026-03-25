CC = mpicxx
CFLAGS = -Wall -O3 --std=c++14
LIBS = -larmadillo 

OBJS = main.o random.o system.o individual.o

EXE = main.exe

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(OBJS) -o $(EXE) $(LIBS)

main.o: main.cpp system.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

individual.o: individual.cpp individual.h random.h
	$(CC) -c individual.cpp -o individual.o $(CFLAGS)

system.o: system.cpp system.h individual.h random.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)

clean:
	rm -f *.o $(EXE) seed.out