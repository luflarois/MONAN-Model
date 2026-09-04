.SUFFIXES: .c .o

all: solv 

OBJS =  spAllocate.o \
		spBuild.o \
		spFactor.o \
		spUtils.o \
		spOutput.o \
		spSolve.o \
		spTest.o \
		spFortran.o

solv: $(OBJS)
	ar -ru lib_solver.a $(OBJS)

# Dependência explícita (opcional, pois o make já infere pelo nome)


# Regra para compilar .c em .o
.c.o:
	$(CC_SERIAL) $(CFLAGS) -c $< $(CINCLUDES) -I./

clean:
	$(RM) *.o *.mod *.i *.a