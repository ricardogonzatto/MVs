all: 
	python srcs/conversor.py
	gcc ./srcs/poli2.c -o ./exec/poli2.o -lm -lgmp -lmps -lpthread
	./exec/poli2.o
