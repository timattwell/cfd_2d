default: all


Model.o: Model.cpp Model.h
	g++ -std=c++11 -Wall -o Model.o -c Model.cpp

Burgers.o: Burgers.cpp Burgers.h Model.h
	g++ -std=c++11 -Wall -o Burgers.o -c Burgers.cpp

compile:  Model.o Burgers.o
	g++ -o my_prog Model.o Burgers.o 

diff: compile
	./my_prog 0.0 0.0 0.0 1.0
	
advx: compile
	./my_prog 1.0 0.0 0.0 0.0

advy: compile
	./my_prog 0.0 1.0 0.0 0.0

burg: compile
	./my_prog 1.0 0.5 1.0 0.02


.PHONY: clean # Specify that ’clean’ is not a real file
	target

clean:
	rm -f *.o my_prog   # Clean up (and ignore any errors)

all: burg clean