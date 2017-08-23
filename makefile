OPTIONS = -O2 -Wall
EXEC = afmCond

all: executable
executable: src_C/afmCond.c
	gcc $(OPTIONS) -o $(EXEC) src_C/afmCond.c -lm 

#verdict: 
#best (-O2) test in c:   34s
#test (-O2) test in c++: 41s
#the end.
