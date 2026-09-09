.SUFFIXES: .f90 .o

# Desabilita as regras implícitas embutidas do GNU Make (ex: regra de
# Modula-2 "%.o: %.mod" via m2c), que colidem com os .mod do Fortran
# e causavam dependência circular em mem_spack.mod <-> mem_spack.o.
MAKEFLAGS += -r

#.PHONY: solver_chem mech tuv_rad

all: solver_chem mech tuv_rad chemistry_monan  

OBJS =  mem_spack.o \
		modChem.o \
		monan_chemistry_vars.o \
		monan_chemistry_interface.o \
		mod_rodas3.o \
		modMemoryChem.o \
		modTest.o \
		mod_chem_output.o

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
monan_chemistry_vars.o: monan_chemistry_vars.f90
modChem.o: modChem.f90 monan_chemistry_vars.o monan_chemistry_interface.o mod_rodas3.o modMemoryChem.o \
	tuv_rad modTest.o mod_chem_output.o
monan_chemistry_interface.o: monan_chemistry_interface.f90 monan_chemistry_vars.o modMemoryChem.o
modMemoryChem.o: modMemoryChem.f90

# Regra para compilar .f90 em .o
.f90.o:
	$(FC) $(FFLAGS) -r -c $< $(FCINCLUDES) -I.. -I../../framework -I../../external/esmf_time_f90 -I./tuv -I./relacs_tuv -I../physics

clean:
	( cd tuv; $(MAKE) clean )
	( cd utils_c; $(MAKE) clean )
	( cd relacs_tuv; $(MAKE) clean )
	$(RM) *.o *.mod *.i *.a
