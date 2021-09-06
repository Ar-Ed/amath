C = g++
opt = -std=c++20 -O3

all: 
	${C} ${opt} lib.cpp test.cpp -o test
	./test
