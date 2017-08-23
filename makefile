OPTIONS = -O2 -Wall
EXEC = afmCond

all: executable
executable: src/afmCond.c
	gcc $(OPTIONS) -o $(EXEC) src/afmCond.c -lm  

