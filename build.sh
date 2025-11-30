#!/bin/sh

g++ ./*.cpp \
	-Wall \
	-O3 \
	-march=native \
	-mavx2 \
	-lgmp \
	-o ./main

#./main
