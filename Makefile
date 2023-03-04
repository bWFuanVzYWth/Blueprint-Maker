test: main.c
	gcc main.c -o test -static -fexec-charset=GBK -Wall -Wfloat-equal -Wshadow -Wbad-function-cast -Wcast-qual -Wcast-align -Wsign-compare -Waggregate-return -Wunreachable-code -Wconversion -Ofast -flto -pipe -march=native -mtune=native -fopt-info -fopenmp

clear:
	rm test.exe