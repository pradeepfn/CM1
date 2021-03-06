
#SHELL = /bin/sh


#-----------------------------------------------------------------------------
#                      NETCDF SECTION
#  (uncomment the following four lines if you want netcdf output capability)
#   (also, make sure the paths to netcdf files are correct for your machine)
#              (NOTE: Don't change lines 3 and 4!)
#
#OUTPUTINC = -I/usr/local/netcdf/include
#OUTPUTLIB = -L/usr/local/netcdf/lib
#OUTPUTOPT = -DNETCDF
#LINKOPTS  = -lnetcdf
#-----------------------------------------------------------------------------
#                         HDF SECTION
#              NEW for cm1r15:  CM1 now supports HDF5
#  (uncomment the following four lines if you want hdf output capability)
#   (also, make sure the paths to hdf files are correct for your machine)
#              (NOTE: Don't change lines 3 and 4!)
# Note: You may need to remove -lsz.
#
#OUTPUTINC = -I/usr/local/hdf5/include
#OUTPUTLIB = -L/usr/local/hdf5/lib -L/usr/local/szip/lib
#OUTPUTOPT = -DHDFOUT
#LINKOPTS  = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lsz -lz -lm
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#                     HARDWARE SECTION
#-- Choose the appropriate architecture, and uncomment all lines 
#-- in that section.  See "README.compile" for more information.
#-----------------------------------------------------------------------------
#  Linux, single processor, Portland Group compiler
#FC   = pgf90
#OPTS = -Mfree -I../include -pc 64 -Kieee -Ktrap=none -O2
#CPP  = cpp -C -P -traditional
#-----------------------------------------------------------------------------
#  Linux, OpenMP shared memory, Portland Group compiler
#FC   = pgf90
#OPTS = -Mfree -I../include -pc 64 -Kieee -Ktrap=none -O2 -mp
#CPP  = cpp -C -P -traditional
#-----------------------------------------------------------------------------
#  Linux, distributed memory, Portland Group compiler (NCAR's lightning)
#FC   = mpif90
#OPTS = -Mfree -I../include -pc 64 -Kieee -Ktrap=none -O2
#CPP  = cpp -C -P -traditional
#DM   = -DMPI
#-----------------------------------------------------------------------------
#  Linux, single processor, Intel fortran compiler (eg, SHARCNET's saw)
#FC   = ifort
#OPTS = -I../include -O3 -FR -assume byterecl -fp-model precise
#CPP  = cpp -C -P -traditional
#-----------------------------------------------------------------------------
#  Linux, distributed memory, Intel compiler (eg, SHARCNET's saw)
#FC   = mpif90
#OPTS = -I../include -O3 -ftz -assume byterecl -fp-model precise
#CPP  = cpp -C -P -traditional
#DM   = -DMPI
#-----------------------------------------------------------------------------
#  Linux, distributed memory, PathScale compiler (eg, SHARCNET's whale)
#FC   = mpif90
#OPTS = -freeform -I../include -O2 -OPT:Olimit=0
#CPP  = cpp -C -P -traditional
#DM   = -DMPI
#-----------------------------------------------------------------------------
#  Linux, single processor, Portland Group compiler, double precision
#FC   = pgf90
#OPTS = -Mfree -I../include -pc 64 -Kieee -Ktrap=none -O2 -r8
#CPP  = cpp -C -P -traditional
#DP   = -DDP
#-----------------------------------------------------------------------------
#  Linux, OpenMP shared memory, Portland Group compiler, double precision
#FC   = pgf90
#OPTS = -Mfree -I../include -pc 64 -Kieee -Ktrap=none -O2 -r8 -mp
#CPP  = cpp -C -P -traditional
#DP   = -DDP
#-----------------------------------------------------------------------------
#  Linux and Mac, using g95 compiler
#FC   = g95
#OPTS = -ffree-form -I../include -O2
#CPP  = cpp -C -P -traditional
#-----------------------------------------------------------------------------
#  Linux and Mac, using g95 compiler, multiple processors, distributed memory using MPI
#OPTS = -f90=g95 -ffree-form -I../include -O2
FC = mpif90 
OPTS = -g -I../include -I/home/pradeep/checkout/tacc_phoenix/include -L/home/pradeep/checkout/tacc_phoenix/lib -lphoenix -larmci
CPP = cpp  -P -traditional
DM = -DMPI
#-----------------------------------------------------------------------------
#  IBM POWER6, single processor
#      (e.g., NCAR's bluefire)
#      NOTE: on bluefire, use 'gmake'
#FC   = xlf_r
#OPTS = -qfree=f90 -I../include -O3 -qstrict -qnozerosize -qmaxmem=-1 -q64
#CPP  = /lib/cpp -C -P
#ADV  = -DSGN
#-----------------------------------------------------------------------------
#  IBM POWER6, multiple processors, shared memory using OpenMP
#      (e.g., NCAR's bluefire)
#      NOTE: on bluefire, use 'gmake'
#FC   = xlf_r
#OPTS = -qfree=f90 -I../include -O3 -qstrict -qnozerosize -qmaxmem=-1 -q64 -qsmp=noauto
#CPP  = /lib/cpp -C -P
#ADV  = -DSGN
#-----------------------------------------------------------------------------
#  IBM POWER6, multiple processors, distributed memory using MPI
#      (e.g., NCAR's bluefire)
#      NOTE: on bluefire, use 'gmake'
#FC   = mpxlf_r
#OPTS = -qfree=f90 -I../include -O3 -qstrict -qnozerosize -qmaxmem=-1 -q64
#CPP  = /lib/cpp -C -P
#DM   = -DMPI
#ADV  = -DSGN
#-----------------------------------------------------------------------------
#  IBM POWER6, multiple processors, hybrid ... both OpenMP and MPI
#      (e.g., NCAR's bluefire)
#      NOTE: on bluefire, use 'gmake'
#FC   = mpxlf_r
#OPTS = -qfree=f90 -I../include -O3 -qstrict -qnozerosize -qmaxmem=-1 -q64 -qsmp=noauto
#CPP  = /lib/cpp -C -P
#DM   = -DMPI
#ADV  = -DSGN
#-----------------------------------------------------------------------------
#  SGI, distributed memory (NCSA's cobalt)
#FC   = ifort
#OPTS = -lmpi -O2 -FR -assume byterecl -ftz -I../include -I/usr/apps/hdf/include
#CPP  = cpp -C -P -traditional
#DM   = -DMPI
#-----------------------------------------------------------------------------
#-- You shouldn't need to change anything below here
#-----------------------------------------------------------------------------

SRC   = allocate.F  adv.F \
	advaxi.F \
	anelp.F \
	base.F \
	bc.F \
	cm1.F \
	comm.F \
	diff2.F \
	diff6.F \
	goddard.F \
	init3d.F \
	init_physics.F \
	init_surface.F \
	init_terrain.F \
	kessler.F \
	lfoice.F \
	maxmin.F \
	misclibs.F \
	morrison.F \
	param.F \
	parcel.F \
	pdef.F \
	poiss.F \
	sfcphys.F \
	singleton.F \
	solve.F \
	sounde.F \
	sound.F \
	soundns.F \
	statpack.F \
        thompson.F \
	turbtke.F \
	writeout.F \
	writeout_nc.F \
	writeout_hdf5.F \
	raddata3d.F \
	radlib3d.F \
	irrad3d.F \
	sorad3d.F \
	radtrns3d.F \
	getcape.F \
	ysu.F \
	sfclay.F \
	slab.F \
	oml.F 

OBJS = $(addsuffix .o, $(basename $(SRC)))

FFLAGS  =  $(OPTS)
AR      = ar cru

.SUFFIXES:
.SUFFIXES:      .F .f90 .o

all : cm1

cm1:			$(OBJS)
			$(FC) $(OBJS) $(FFLAGS) $(OUTPUTINC) $(OUTPUTLIB) $(LINKOPTS) -o ../run/cm1.exe

.F.o:
			$(CPP) $(DM) $(DP) $(ADV) $(OUTPUTOPT) $*.F > $*.f90
			$(FC) $(FFLAGS) $(OUTPUTINC) -c $*.f90

code:
			$(AR) onefile.F $(SRC)
			mv onefile.F ../run

clean:
			rm -f *.f90 *.o *.a *.mod

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

poiss.o: singleton.o
param.o: thompson.o
param.o: morrison.o
solve.o: thompson.o morrison.o sfclay.o ysu.o oml.o slab.o
init_physics.o: sfclay.o oml.o slab.o
init_surface.o: oml.o
