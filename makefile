C = g++
opt = -std=c++20 

all: 
	${C} ${opt} array.cpp test.cpp -o test -I /usr/local/Cellar/python@3.9/3.9.6/Frameworks/Python.framework/Versions/3.9/include/python3.9 \
 	-I /Library/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages/numpy/core/include \
 	-L /usr/local/Cellar/python@3.9/3.9.6/Frameworks/Python.framework/Versions/3.9/lib \
	-lpython3.9
	./test

cmp:
	${C} ${opt} complex.cpp testcomp.cpp -o testc
	./testc
