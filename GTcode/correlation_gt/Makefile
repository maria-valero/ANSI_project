CC = gcc
CFLAGS = -g
CFLAGS = -pg
CFLAGS = -O3

low_pass_downsampling: low_pass_downsampling.o filtersubs.o
	gcc ${CFLAGS} -o low_pass_downsampling low_pass_downsampling.o filtersubs.o -lpthread -lm
