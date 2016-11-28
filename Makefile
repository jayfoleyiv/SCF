LIB        = -L. -lm -llapack -lblas -Lg2c
INCLUDE    = -I.
CFLAGS     = -O2 -g -pg
EXEC       = HF.x
CXX        = g++

${EXEC}: HF.c  blas.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} HF.c blas.o -o ${EXEC}


blas.o: blas.c blas.h
	${CXX} ${LIB} -c blas.c ${CFLAGS}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<

