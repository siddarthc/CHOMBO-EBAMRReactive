# compile and link makefile for CHEMKIN-II related  transport codes
# JES Oct 2003
#

FF = ifort 
FOPT = -O2

# Transport property evaluation


dtp.o :	dtp.f                   
	$(FF) $(FOPT) -c dtp.f

dtp: dtp.o  cklib.o dmath.o  
	$(FF) -o dtp dtp.o cklib.o dmath.o tranlib.o

# use flag below for building dmath.o in case compiling aborts due to 
# precission error (see [info -f g77 M GLOBALS])
# fixed dmath.o is fixed, original file in dmath.o.original
# pintgen 28.10.03
# dmath.o: dmath.f
#	$(FF) -fno-globals -c dmath.f

dmath.o: dmath.f
	$(FF) $(FOPT) -c dmath.f

ckinterp: ckinterp.o
	$(FF) -o ckinterp ckinterp.o

ckinterp.o: ckinterp.f
	$(FF) $(FOPT) -c ckinterp.f

cklib.o: cklib.f
	$(FF) $(FOPT) -c cklib.f

tranlib.o: tranlib.f
	$(FF) $(FOPT) -c tranlib.f

tranfit.o: tranfit.f
	$(FF) $(FOPT) -c tranfit.f

tranfit: tranfit.o dmath.o tranlib.o cklib.o
	$(FF)  -o tranfit tranfit.o dmath.o tranlib.o cklib.o

# Thermal conductivity and Lewis number evaluation

tcon.o :	tcon.f                    # main program
	$(FF) $(FOPT) -c tcon.f

tcon:   tcon.o  cklib.o dmath.o
	$(FF) -o tcon tcon.o cklib.o dmath.o tranlib.o

# simple version of transport evaluation for mixture 

tc.o :	tc.f                    # main program
	$(FF) $(FOPT) -c tc.f

tc : tc.o   cklib.o  dmath.o
	$(FF) -o tc tc.o cklib.o dmath.o tranlib.o

#####################################################################
# helpers to build everything from scratch (*.f's)

clean:
	rm -f tcon tcon.o tc tc.o dtp dtp.o tran.bin tranfit \
	dmath.o tranfit.o tranlib.o cklib.o ckinterp ckinterp.o \
	chem.bin chem.out tran.inp tcon.out

all:    ckinterp.f ckstrt.f tranlib.f chem.inp dmath.f dtp.f tcon.f \
	therm.dat tranfit.f cklib.f tc.f tran.dat
	make ckinterp.o; make ckinterp; ./ckinterp; \
	make cklib.o; make tranlib.o; make dmath.o; make tranfit.o; \
	make tranfit; ./tranfit; make dtp.o; \
	make dtp; make tcon.o; make tcon; make tc.o; make tc;  
