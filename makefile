C = g++
opt = -std=c++20 -O3

all: 
	${C} ${opt} array.cpp test.cpp -o test
	./test

cmp:
	${C} ${opt} complex.cpp testcomp.cpp -o testc
	./testc

ra:
	${C} ${opt} rational.cpp testrat.cpp -o testr
	./testr