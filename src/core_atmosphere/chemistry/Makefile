.SUFFIXES: .f90 .o

all: solver_chem mech tuv_rad chemistry_monan  

OBJS =  mem_spack.o \
		modChem.o \
		monan_chemistry_vars.o \
		monan_chemistry_interface.o \
		mod_rodas3.o \
		modMemoryChem.o

chemistry_monan: $(OBJS)
	ar -ru lib_chem.a $(OBJS)

tuv_rad: 
	( cd tuv; $(MAKE) all) 
	( echo "****** compiling tuv_rad ******" )

solver_chem:
	( cd utils_c; $(MAKE) all) 
	( echo "****** compiling c utils solver ******" )

mech:
	( cd relacs_tuv; $(MAKE) all) 
	( echo "****** compiling mech relacs_tuv******" )

# Dependência explícita (opcional, pois o make já infere pelo nome)
mod_rodas3.o: mod_rodas3.f90 modMemoryChem.o
modChem.o: monan_chemistry_vars.o monan_chemistry_interface.o mod_rodas3.o modMemoryChem.o
monan_chemistry_interface.o: monan_chemistry_interface.f90 monan_chemistry_vars.o

# Regra para compilar .f90 em .o
.f90.o:
	$(FC) $(FFLAGS) -r -c $< $(FCINCLUDES) -I.. -I../../framework -I../../external/esmf_time_f90 -I./tuv -I./relacs_tuv -I../physics

clean:
	$(RM) *.o *.mod *.i