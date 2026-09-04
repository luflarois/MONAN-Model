.SUFFIXES: .f90 .o

all: relacs_tuv  

OBJS =  chem_list.o \
		chem_spack_dratedc.o \
		chem_spack_fexchem.o \
		chem_spack_fexloss.o \
		chem_spack_fexprod.o \
		chem_spack_jacdchemdc.o \
		chem_spack_kinetic.o \
		chem_spack_rates.o 

relacs_tuv: $(OBJS)
	ar -ru lib_relacs_tuv.a $(OBJS)

# Dependência explícita (opcional, pois o make já infere pelo nome)
chem_spack_fexchem.o: chem_spack_rates.o

# Regra para compilar .f90 em .o
.f90.o:
	$(FC) $(FFLAGS) -c $< $(FCINCLUDES) -I.. 

clean:
	$(RM) *.o *.mod *.i *.a