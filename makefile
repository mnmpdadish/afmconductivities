OPTIONS = -O2 -Wall -ffast-math
EXEC = afmCond

all: executable hall
executable:
	gcc $(OPTIONS) -o $(EXEC) src/afmCond.c -lm  

hall: 
	gcc $(OPTIONS) -o afmHall src/afmHall.c -lm  

