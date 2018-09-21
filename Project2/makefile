# Comment lines
# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= c++ -Wall
PROG= tridiagMatrix
# this is the math library in C, not necessary for C++
LIB = -lm -llapack -lblas -larmadillo
# Here we make the executable file
${PROG} :	${PROG}.o
	${CC} ${PROG}.o ${LIB} -o ${PROG}
# whereas here we create the object file
${PROG}.o :	${PROG}.cpp
	${CC} -c ${PROG}.cpp
