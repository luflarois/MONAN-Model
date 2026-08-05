.SUFFIXES: .f90 .o

all: tuv_rad

OBJS =  ModTuv.o  \
		ModTuv_driver.o 


tuv_rad: $(OBJS)
	ar -ru lib_tuv.a $(OBJS)

# Dependência explícita (opcional, pois o make já infere pelo nome)
ModTuv.o: ModTuv.f90 #chem_list.o
ModTuv_driver.o: ModTuv_driver.f90 ModTuv.o #chem_list.o

# Regra para compilar .f90 em .o
.f90.o:
	$(FC) $(FFLAGS) -c $< $(FCINCLUDES) -I.. -I../relacs_tuv -I../../../framework -I../../external/esmf_time_f90 -I../../physics -I../

clean:
	$(RM) *.o *.mod *.i