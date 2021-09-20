C = g++
opt = -std=c++20 -O0


all:
	${C} ${opt} array.cpp test.cpp -I/Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/include/python2.7 -lpython2.7 -pthread -o test
	./test

rm:
	rm test
	