# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(RAD_DIR)_$(MICRO_DIR)

# Determine platform 
PLATFORM := $(shell uname -s)

#----------------------------------
# AIX (tested only on IBM SP)
#

ifeq ($(PLATFORM),AIX)

#INC_MPI      := /usr/local/include
#LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib


FF77 = mpxlf90_r -c -qsuffix=f=f -qfixed=132
FF90 = mpxlf90_r -c -qsuffix=f=f90
CC = cc -c -DAIX
FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q -I${INC_NETCDF}
#FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en -I${INC_NETCDF}
LD = mpxlf90_r
LDFLAGS = -bmaxdata:512000000 -bmaxstack:256000000 -L${LIB_NETCDF} -lnetcdf

endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(PLATFORM),IRIX64)

INC_MPI      := /usr/local/include
LIB_MPI      := /usr/local/lib
INC_NETCDF   := /usr/local/include
LIB_NETCDF   := /usr/local/lib

FF77 = f90 -c -fixedform  -extend_source
FF90 = f90 -c -freeform
CC = cc -c -DIRIX64
FFLAGS = -O3 
#FFLAGS = -g -DEBUG:subscript_check=ON:trap_uninitialized=ON 
FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
LD = f90 
LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

endif
#----------------------------------------------------------------------
# Linux, Intel Compiler
#

#ifeq ($(PLATFORM),Linux)
#
#LIB_MPI = /usr/local/pkg/iopenmpi/lib
#INC_MPI = /usr/local/pkg/iopenmpi/include
#INC_NETCDF = /nfs/user08/marat/local/include
#LIB_NETCDF = /nfs/user08/marat/local/lib
#
#
#FF77 = /usr/local/pkg/iopenmpi/bin/mpif90 -c -fixed -extend_source
#FF90 = /usr/local/pkg/iopenmpi/bin/mpif90 -c
#CC = mpicc -c -DLINUX
#
#
#FFLAGS = -O3 
##FFLAGS = -g -ftrapuv -check all
#
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
#LD = /usr/local/pkg/iopenmpi/bin/mpif90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf
#

#endif
#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at San Diego SC
#
#ifeq ($(PLATFORM),Linux)

#INC_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/include
#LIB_NETCDF   := /usr/local/apps/V1R3/netcdf-3.6.0-p1/lib

#FF77 = mpxlf90  -qarch=440 -qsuffix=f=f -qfixed=132
#FF90 = mpxlf90  -qarch=440 -qsuffix=f=f90
#CC = mpcc -c -DLinux
#FFLAGS = -c -O3 -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF}
#LD = mpxlf90
#LDFLAGS = -L${LIB_NETCDF} -lnetcdf

#endif
#----------------------------------------------------------------------
# Linux, XLF compiler, Bluegene at BNL
#
#
#ifeq ($(PLATFORM),Linux)
#
#ifeq ($(HOST),fen)
## Compute nodes:
#INC_NETCDF   := /gpfs/home2/marat/local/netcdf/include
#LIB_NETCDF   := /gpfs/home2/marat/local/netcdf/lib
#FF77 = mpixlf90 -qarch=440 -qsuffix=f=f -qfixed=132
#FF90 = mpixlf90 -qarch=440 -qsuffix=f=f90
#CC = mpicc -c -DLinux
#FFLAGS = -c -O3 -g -qtune=440 -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF}
#LD = mpixlf90  
#LDFLAGS = -L${LIB_NETCDF} -L/bgl/local/lib/ -L/opt/ibmmath/lib -lnetcdf -llapack.rts -lesslbg
#
#else
## front end nodes and vis nodes
#INC_NETCDF   := /gpfs/home2/marat/netcdf-vis-3.6.2/include
#LIB_NETCDF   := /gpfs/home2/marat/netcdf-vis-3.6.2/lib
#LIB_MPI = /usr/local/mpich2-1.0.5p4/lib
#INC_MPI = /usr/local/mpich2-1.0.5p4/include
#FF77 = xlf90 -qsuffix=f=f -qfixed=132
#FF90 = xlf90 -qsuffix=f=f90
#CC = cc -c -DLinux
#FFLAGS = -c -O3  -qstrict -qmaxmem=-1 -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS +=  -I${INC_NETCDF} -I${INC_MPI}
#LD = xlf90 
#LDFLAGS = -L${LIB_NETCDF} -L${LIB_MPI} -L/gpfs/home1/slatest/blas/BLAS  -L/gpfs/home1/slatest/vislapack-2/lapack-3.2 -lmpich -lnetcdf -llapack_RS6K -lblas_LINUX 
#endif
#
#
#endif
#----------------------------------------------------------------------
# Linux, Intel Compiler, habaneros cluster at Columbia
#

# AAW, modified 11/16/2016

ifeq ($(PLATFORM),Linux)
LIB_MPI = /rigel/opt/parallel_studio_xe_2017/impi/2017.0.098/intel64/lib
INC_MPI = /rigel/opt/parallel_studio_xe_2017/impi/2017.0.098/intel64/include
LIB_NETCDF = /rigel/opt/netcdf-fortran-4.4.4/lib
INC_NETCDF = /rigel/opt/netcdf-fortran-4.4.4/include
#
#
#FF77 = mpif90 -c -fixed -extend_source
#FF90 = mpif90 -c
#CC = mpicc -c -DLINUX
#FF77 = mpiifort -c -fixed_extend_source
FF77 = mpiifort -c -fixed -extend_source
FF90 = mpiifort -c
CC = mpiicc -c -DLINUX
#
#
FFLAGS = -O3
##FFLAGS = -g -ftrapuv -check all
#
#FFLAGS += -I${INC_MPI} -I${INC_NETCDF} #normally this is uncommented
FFLAGS += -g -I${INC_MPI} -I${INC_NETCDF}
#LD = mpif90
LD = mpiifort
#LDFLAGS = -L${LIB_NETCDF} -lnetcdff -lnetcdf -lnetcdf
LDFLAGS = -Wl,-rpath=${LIB_NETCDF},--enable-new-dtags -L/rigel/opt/parallel_studio_xe_2017/impi/2017.0.098/intel64/lib -L${LIB_NETCDF} -L${LIB_MPI} -lifcore -lmpi -lnetcdff
#

#endif
#----------------------------------------------------------------------
# Linux, Portland Group Compiler
#
#----------------------------------------------------------------------
# Linux, Portland Group Compiler
#

# bloss:
#

#
# -- twc 9/2011 -- modified for use on elorenz
#  adapted (/copied) from work by Alli
#
#
#ifeq ($(PLATFORM),Linux)
#
#  # Default compiler flags
#  FFLAGS =  -Kieee -fastsse #-g -C -Ktrap=fp #
#
#  ifeq ($(HOSTNAME),olympus)
#    # special setup on olympus
#    MPIF90 = /usr/local/openmpi/bin/mpif90
#    FFLAGS += -tp k8-64
#
#  else
#    # Determine platform
#    PROCESSOR := $(shell uname -p)
#
#    # DEFAULT MPI LOCATION AND COMPILER
##    LIB_MPI = /usr/local/mpich/lib
##    INC_MPI = /usr/local/mpich/include
#    MPIF90 = mpif90

#    # CHANGE COMPILER IF RUNNING ON 64BIT MACHINE.
#    ifeq ($(PROCESSOR),x86_64)
##      MPIF90 = /usr/local/pgi/linux86-64/7.1-1/bin/pgf90
#      FFLAGS += -tp k8-64
#    else
##      # ADD -tp k8-32 IF RUNNING ON REX.
#      ifeq ($(HOSTNAME),rex)
#        FFLAGS += -tp k8-32
#      else
#        FFLAGS += -tp k7
#      endif
#    endif
#
#    ifeq ($(HOSTNAME),rex)
#      # UNCOMMENT TO USE LAHEY COMPILER -- USEFUL FOR DEBUGGING
#      LIB_MPI = /usr/local/mpich-lf/lib
#      INC_MPI = /usr/local/mpich-lf/include
#      MPIF90 = lf95
#      FFLAGS = -g --trap --chk aesu #--trap # --o2 #
#    endif
#
##    # UNCOMMENT FOR MYRINET RUN
##    #LIB_MPI = /usr/local/mpich-gm/lib
##    #INC_MPI = /usr/local/mpich-gm/include
##    #MPIF90 = /usr/local/mpich-gm/bin/mpif90
##
##    # ADD MPI FLAGS HERE IF NOT USING mpif90
##    FFLAGS += -I${INC_MPI}
##    LDFLAGS += -L${LIB_MPI} -lmpich
#  endif
#
#  FF77 = ${MPIF90} -c
#  FF90 = ${MPIF90} -c
#  CC = gcc -c -DLINUX -g
#
#  LD = ${MPIF90} ${FFLAGS}
#  FFLAGS += -I${INC_NETCDF}
#  LDFLAGS += -L${LIB_NETCDF} -L{LIB_MPI} -lnetcdf -lmpi -lnetcdff

## end bloss:

# older options:

#LIB_MPI = /usr/pgi/linux86/5.1/lib
#INC_MPI = /usr/pgi/linux86/5.1/include
#LIB_MPI = /usr/local/lam-mpi/lib
#INC_MPI = /usr/local/lam-mpi/include

#FF77 = pgf90 -c
#FF90 = pgf90 -c -Mfreeform
#CC = cc -c  -DLINUX

#FFLAGS = -Mnoframe -Mvect -Munroll -O2 -Mbyteswapio  

#FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
#LD = pgf90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf

endif

#--------------------------------------------
# Apple Mac OS X (Darwin) (Absoft Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/absoft/include
#LIB_NETCDF   := /usr/local/absoft/lib
#INC_MPI      := /usr/local/absoft/include
#LIB_MPI       := /usr/local/absoft/lib

#FF77 = f90 -c -f fixed
#FF90 = f90 -c -f free
#CC = cc -c -DMACOSX

#FFLAGS = -O3 -noconsole -nowdir -YEXT_NAMES=LCS -s -YEXT_SFX=_  -z4
#LD = f90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf

#endif

#--------------------------------------------
# Apple Mac OS X (Darwin) (NAG Fortran)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/nag/include
#LIB_NETCDF   := /usr/local/nag/lib
#INC_MPI      := /usr/local/nag/include
#LIB_MPI       := /usr/local/nag/lib

#FF77 = f95 -c -fixed -kind=byte
#FF90 = f95 -c -free -kind=byte
#CC = cc -c -DMACOSX

#FFLAGS =      # don't use any optimization -O* option! Will crash!
##FFLAGS = -gline -C=all -C=undefined  # use for debugging
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = f95 
##LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpich -lnetcdf 
#LDFLAGS =  

#endif


#----------------------------------
# Apple Mac OS X (Darwin) (XLF compiler)
#

#ifeq ($(PLATFORM),Darwin)

#INC_NETCDF   := /usr/local/xlf/include
#LIB_NETCDF   := /usr/local/xlf/lib
#INC_MPI      := /usr/local/xlf/include
#LIB_MPI       := /usr/local/xlf/lib

#FF77 = xlf90 -c -qsuffix=f=f -qfixed=132
#FF90 = xlf90 -c -qsuffix=f=f90
#CC = cc -c -DMACOSX 
#FFLAGS = -c -O3 -qstrict -qmaxmem=-1 -qarch=auto -qspillsize=5000 -Q
##FFLAGS = -c -qinitauto=FF -g -qflttrap=zerodivide:enable -qflttrap=ov:zero:inv:en
#FFLAGS += -I$(SAM_SRC)/$(RAD_DIR) -I${INC_MPI} -I${INC_NETCDF}
#LD = xlf90
#LDFLAGS = -L${LIB_MPI} -L${LIB_NETCDF} -lmpi -lnetcdf

#endif

#----------------------------------

#----------------------------------
# Apple Mac OS X (Darwin) (Intel compiler)
#

ifeq ($(PLATFORM),Darwin)

INC_MPI      := /usr/local/include
LIB_MPI       := /usr/local/lib
INC_NETCDF      := /usr/local/include
LIB_NETCDF       := /usr/local/lib


FF77 = mpif90 -c -fixed -extend_source
FF90 = mpif90 -c 
CC = mpicc -c -DLINUX


FFLAGS = -O3 -pad
#FFLAGS = -g -ftrapuv -check all

FFLAGS += -I${INC_MPI} -I${INC_NETCDF}
LD = mpif90 
LDFLAGS = -L${LIB_NETCDF} -lnetcdf 

endif


#----------------------------------
#----------------------------------------------
# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<



include Depends



clean: 
	rm ./OBJ/*


