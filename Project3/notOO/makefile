#compiler options and target
CC= c++ -Wall
PROG= earthsunsystem
#libraries
LIB = -lm -llapack -lblas -larmadillo
#executable file
${PROG}.exe :	${PROG}.o
	${CC} ${PROG}.o ${LIB} -o ${PROG}.exe
#object file
${PROG}.o :	${PROG}.cpp
	${CC} -c ${PROG}.cpp

