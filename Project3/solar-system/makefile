
#compiler options and target
CXX= c++ -Wall
PROG= main

#libraries
LIB = -lm -llapack -lblas -larmadillo 
#HEADERS= solarsystem.h celestialbody.h vec3.h

OBJECTS= solarsystem.o celestialbody.o vec3.o solver.o

#default goal:
${PROG}.exe :	${PROG}.o ${OBJECTS}
	${CXX} ${PROG}.o ${LIB} -o ${PROG}.exe ${OBJECTS}


#other recipes
${PROG}.o :	${PROG}.cpp ${PROG}.h
	${CXX} -c ${PROG}.cpp

solarsystem.o : solarsystem.cpp solarsystem.h
	${CXX} -c solarsystem.cpp

celestialbody.o : celestialbody.cpp celestialbody.h
	${CXX} -c celestialbody.cpp

vec3.o : vec3.cpp vec3.h
	${CXX} -c vec3.cpp

solver.o : solver.cpp solver.h
	${CXX} -c solver.cpp
