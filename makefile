OPTIONS = -O2 -Wall -ffast-math
EXEC = afmCond

all: executable
executable: src/afmCond.c
	gcc $(OPTIONS) -o $(EXEC) src/afmCond.c -lm  

