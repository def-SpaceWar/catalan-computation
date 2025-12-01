#!/bin/sh

g++ ./*.cpp \
	-Wall \
	-Ofast \
	-march=native \
	-mavx2 \
	-lgmp \
	-o ./main

#./main
