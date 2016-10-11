LIB        = -L. -lm -llapack -lblas -Lg2c
INCLUDE    = -I.
CFLAGS     = -O2 -g -pg
EXEC       = main.x
CXX        = g++

${EXEC}: main.c memory.o blas.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} main.c blas.o memory.o -o ${EXEC}


memory.o: memory.c memory.h
	${CXX} -c memory.c ${CFLAGS}

blas.o: blas.c blas.h
	${CXX} ${LIB} -c blas.c ${CFLAGS}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

