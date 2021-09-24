C = g++
opt0 = -std=c++20 -O0
opt3 = -std=c++20 -O3


all:
	${C} ${opt0} test.cpp -I/Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/include/python2.7 -lpython2.7 -pthread -o test
	./test

pre:##??
	${C} ${opt0} -x c++-header ./HEADER_ONLY_AMATH/amath.h -I/Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/include/python2.7 -lpython2.7 -pthread  

rm:
	rm test
	