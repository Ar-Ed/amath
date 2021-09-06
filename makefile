C = g++
opt = -std=c++20 -O3

all: 
	${C} ${opt} array.cpp test.cpp -o test
	./test
