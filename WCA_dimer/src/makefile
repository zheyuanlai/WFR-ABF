#------------------
#   Compiler
#------------------

CC = g++

#------------------
#   Libraries
#------------------

ITINCLUDE3  =   #-I/... # when include files are requested
ITINCLUDE   = ${ITINCLUDE3}
ITLIB       = -L/usr/lib/ # directory where the required libraries are

#----------------------------
#   Compilation options
#----------------------------

INCLUDE_PATH = -I./ ${ITINCLUDE}
LDFLAGS      = ${ITLIB}
LIBRARY      = -lstdc++ -lc -lm 
WARN         = -Wno-deprecated 
CPPFLAGS     = -ffast-math -funroll-loops ${WARN} ${INCLUDE_PATH}
DEFINES      = -Wall

.cpp.o: ; ${CC} ${DEFINES} ${INCLUDE_PATH} ${WARN} -g -c $*.cpp

#----------------------------
#   Files to be compiled
#----------------------------

OBJsa = gaussian.o input_sampling.o matrix.o vector.o hamiltonian.o sampling.o main.o
OBJmb = gaussian.o input_mbar.o     matrix.o vector.o hamiltonian.o mbar.o     main.o 
OBJab = gaussian.o input_ABF.o      matrix.o vector.o hamiltonian.o ABF.o      main.o
OBJti = gaussian.o input_TI.o       matrix.o vector.o hamiltonian.o TI.o       main.o
OBJjz = gaussian.o input_jarz.o     matrix.o vector.o hamiltonian.o jarz.o     main.o


sampling:	${OBJsa}
		${CC} ${CPPFLAGS} ${DEFINES} ${OBJsa}  -o sampling++ ${LDFLAGS} ${LIBRARY}

mbar:		${OBJmb}
		${CC} ${CPPFLAGS} ${DEFINES} ${OBJmb}  -o mbar++ ${LDFLAGS} ${LIBRARY}

abf:		${OBJab}
		${CC} ${CPPFLAGS} ${DEFINES} ${OBJab}  -o abf++ ${LDFLAGS} ${LIBRARY}

ti:		${OBJti}
		${CC} ${CPPFLAGS} ${DEFINES} ${OBJti}  -o ti++ ${LDFLAGS} ${LIBRARY}

jarz:		${OBJjz}
		${CC} ${CPPFLAGS} ${DEFINES} ${OBJjz} -o jarz++ ${LDFLAGS} ${LIBRARY}


clean:
	/bin/rm -rf *.o *.exe *~ *++ ../output/* ../input/*~


