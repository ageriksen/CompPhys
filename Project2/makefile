# makefile for JacobiMethodEigenvalues
#
# compiler option, libraries and target
CC = c++ -Wall
PROG = JacobiMethodEigenvalues
#Libraries to include:
LIB = -lm -llapack -lblas -larmadillo
#executable
${PROG}.exe :	${PROG}.o
	${CC} ${PROG}.o ${LIB} -o ${PROG}.exe
#Object file
${PROG}.o : ${PROG}.cpp
	${CC} -c ${PROG}.cpp
