########## User-definable stuff ##########
#
###Compiler and compilation options
COMP = gcc
OPTIONS = -Wall -O3 -std=c99 -g
#OPTIONS += -D_HAS_OMP -fopenmp
#OPTIONS += -D_DEBUG
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
GSL_INC = -I/home/damonge/include
GSL_LIB = -L/home/damonge/lib
#COSMOMAD
CSM_INC = -I/home/anze/local/include
CSM_LIB = -L/home/anze/local/lib
#
FFTW_INC = 
FFTW_LIB = -lfftw3
########## End of user-definable ##########

OPTIONS += $(DEFINEFLAGS)

INC_ALL = -I./src $(GSL_INC) $(CSM_INC) $(DAM_INC) $(FFTW_INC)
LIB_ALL = $(GSL_LIB) $(CSM_LIB) -lcosmomad -lgsl -lgslcblas -lm $(FFTW_LIB)

COMMONO = src/common.o
IOO = src/io.o
COSMOO = src/cosmo.o
TRANSFERO = src/transfers.o
SPECTRAO = src/spectra.o
LOGNORMO = src/lognorm.o
FFTLOGO = src/fftlog.o
MAIN = src/limberjack.c
OFILES = $(COMMONO) $(IOO) $(COSMOO) $(TRANSFERO) $(SPECTRAO) $(LOGNORMO) $(FFTLOGO)

EXE = LimberJack

default : $(EXE)

%.o : %.c
	$(COMP) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(EXE) : $(OFILES)
	$(COMP) $(OPTIONS) $(INC_ALL) $(OFILES) $(MAIN) -o $(EXE) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner : 
	rm -f *~ src/*.o src/*~ $(EXE)
