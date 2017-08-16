HEADERS = src/utilities.h
OPTIONS = -O2 -fpermissive -Wall -std=c++11 -DBOUND_CHECK
EXEC = afmCond

all: executable
executable: src/afmCond.cpp
	gcc $(OPTIONS) -o $(EXEC) src/afmCond.cpp -llapack -lblas -lstdc++

