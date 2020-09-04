#******************************************************************************#
#*                       Makefile for FD3Dspher package                       *#
#*                                                                            *#
#*            Author: Wei ZHANG     Email: zhangwei.zw@gmail.com              *#
#*                       Copyright (C) 2008 Wei ZHANG                         *#
#******************************************************************************#

# $Date: 2009-01-16 13:03:38 -0500 (Fri, 16 Jan 2009) $
# $Revision: 73 $
# $LastChangedBy: zhangw $

#######################################################################
#                 Compiler and macro flags setting                    #
#######################################################################

#WHEREAMI := nw8440
WHEREAMI := uri
#WHEREAMI := pku
#WHEREAMI := ustc
#WHEREAMI := generic

# Please choose COMPILER in the included Makefile.opt.$(WHEREAMI)

#DEBUG := ON
STATIC := ON
#GETARG := ON
#VERBOSE :=ON
#CheckOverFlow := ON

#MPITOPO1D :=ON
#MPIBARRIER := ON
#USEOMP :=ON
#MPIBuffered := ON

USEPML := ON
#WITHQS := ON
#DataTypeDouble := ON

SrcSmooth :=ON
SrcSurface := ON

#CondFreeCharac := ON
CondFreeTIMG := ON
#CondFreeVHOC := ON
CondFreeVLOW := ON
#CondFreeVEXT := ON

DFLAG_LIST := DEBUG STATIC GETARG VERBOSE CheckOverFlow \
			  MPITOPO1D USEOMP MPIBARRIER MPIBuffered DataTypeDouble \
              USEPML WITHQS SrcSmooth SrcSurface \
			  CondFreeCharac CondFreeTIMG \
              CondFreeVHOC CondFreeVEXT CondFreeVLOW

DFLAGS := $(foreach flag,$(DFLAG_LIST),$(if $($(flag)),-D$(flag),)) $(DFLAGS)
DFLAGS := $(strip $(DFLAGS))

#######################################################################
#            directories, source files and target names               #
#######################################################################
skeldirs := OBJ bin src
OBJDIR := ./OBJ
BINDIR := ./bin
SRCDIR := ./src
FPPDIR := ./srcF

SRC_MOD = mod_constants.f90 mod_math.f90 mod_string.f90 mod_nfseis.f90 \
		  mod_para.f90 mod_mpi.f90  \
		  mod_grid.f90 mod_io.f90 mod_media.f90 mod_macdrp.f90 mod_src.f90 \
		  mod_custom.f90
OBJ_MOD := $(foreach file,$(SRC_MOD),$(OBJDIR)/$(file:.f90=.o))

SRC_GRID    := seis3d_grid.f90
EXE_GRID    := seis3d_grid
OBJ_GRID    := $(foreach file,$(SRC_GRID),$(OBJDIR)/$(file:.f90=.o))

SRC_METRIC     := seis3d_metric.f90
EXE_METRIC     := seis3d_metric
OBJ_METRIC     := $(foreach file,$(SRC_METRIC),$(OBJDIR)/$(file:.f90=.o))

SRC_MEDIA   := seis3d_media.f90
EXE_MEDIA   := seis3d_media
OBJ_MEDIA   := $(foreach file,$(SRC_MEDIA),$(OBJDIR)/$(file:.f90=.o))
SRC_MEDIA_MPI := seis3d_media_mpi.f90
EXE_MEDIA_MPI := seis3d_media_mpi
OBJ_MEDIA_MPI := $(foreach file,$(SRC_MEDIA_MPI),$(OBJDIR)/$(file:.f90=.o))

SRC_SOURCE  := seis3d_source.f90
EXE_SOURCE  := seis3d_source
OBJ_SOURCE  := $(foreach file,$(SRC_SOURCE),$(OBJDIR)/$(file:.f90=.o))

SRC_STATION := seis3d_station.f90
EXE_STATION := seis3d_station
OBJ_STATION := $(foreach file,$(SRC_STATION),$(OBJDIR)/$(file:.f90=.o))

SRC_WAVE := $(if $(USEPML),mod_abs_npml.f90,mod_abs_exp.f90) \
          seis3d_wave.f90
EXE_WAVE    := seis3d_wave_mpi
OBJ_WAVE    := $(foreach file,$(SRC_WAVE),$(OBJDIR)/$(file:.f90=.o))

SRC_EXPTSNAP   := tool_expt_snap.f90
EXE_EXPTSNAP   := tool_expt_snap
OBJ_EXPTSNAP   :=  $(foreach file,$(SRC_EXPTSNAP),$(OBJDIR)/$(file:.f90=.o))

SRC_EXPTSEISMO := tool_expt_seismo.f90
EXE_EXPTSEISMO := tool_expt_seismo
OBJ_EXPTSEISMO :=  $(foreach file,$(SRC_EXPTSEISMO),$(OBJDIR)/$(file:.f90=.o))

vpath %.F90 $(FPPDIR)

FPP := /usr/bin/cpp
FPPFLAGS := -P -traditional $(foreach dir,$(FPPDIR),-I$(dir)) $(DFLAGS) 

#######################################################################
#                     compiler and option                             #
#######################################################################
include Makefile.opt.$(WHEREAMI)

#######################################################################
#                            target                                   #
#######################################################################
phony_list := skel all preprocess solver postprocess
.PHONY: $(phony_list)

all: skel preprocess solver kernel
#all: skel preprocess solver postprocess kernel

preprocess: $(BINDIR)/$(EXE_GRID)  \
            $(BINDIR)/$(EXE_METRIC)  \
            $(BINDIR)/$(EXE_MEDIA) $(BINDIR)/$(EXE_MEDIA_MPI) \
            $(BINDIR)/$(EXE_SOURCE) $(BINDIR)/$(EXE_STATION)
solver:  $(BINDIR)/$(EXE_WAVE)
#postprocess: $(BINDIR)/$(EXE_EXPTSNAP) $(BINDIR)/$(EXE_EXPTSEISMO)
skel:
	@mkdir -p $(skeldirs)

$(BINDIR)/$(EXE_WAVE): $(OBJ_MOD)  $(OBJ_WAVE)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_WAVE) $(LDFLAGS)
$(BINDIR)/$(EXE_GRID): $(OBJ_MOD) $(OBJ_GRID)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_GRID) $(LDFLAGS)
$(BINDIR)/$(EXE_METRIC): $(OBJ_MOD) $(OBJ_METRIC)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_METRIC) $(LDFLAGS)
$(BINDIR)/$(EXE_MEDIA): $(OBJ_MOD) $(OBJ_MEDIA)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_MEDIA) $(LDFLAGS)
$(BINDIR)/$(EXE_MEDIA_MPI): $(OBJ_MOD) $(OBJ_MEDIA_MPI)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_MEDIA_MPI) $(LDFLAGS)
$(BINDIR)/$(EXE_SOURCE): $(OBJ_MOD) $(OBJ_SOURCE)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_SOURCE) $(LDFLAGS)
$(BINDIR)/$(EXE_STATION): $(OBJ_MOD) $(OBJ_STATION)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_STATION) $(LDFLAGS)
$(BINDIR)/$(EXE_EXPTSNAP): $(OBJ_MOD) $(OBJ_EXPTSNAP)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_EXPTSNAP) $(LDFLAGS)
$(BINDIR)/$(EXE_EXPTSEISMO): $(OBJ_MOD) $(OBJ_EXPTSEISMO)
	$(FC) -o $@ $(OBJ_MOD) $(OBJ_EXPTSEISMO) $(LDFLAGS)

RM := rm
cleanexe:
	$(RM) -f $(BINDIR)/*
cleanobj:
	$(RM) -f $(OBJDIR)/*
cleanf90:
	$(RM) -f $(SRCDIR)/*
cleanall: cleanexe cleanobj cleanf90
distclean: cleanexe cleanobj cleanf90

#######################################################################
#                        subffixes rules                              #
#######################################################################

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

seis3d_media_mpi.f90: seis3d_media.F90
	$(FPP) $(FPPFLAGS) -DMediaMPI $< > $(SRCDIR)/$(@F)

seis3d_wave.f90: seis3d_wave_8op_p3.F90
	$(FPP) $(FPPFLAGS) $< > $(SRCDIR)/$(@F)

%.f90 : %.F90
	$(FPP) $(FPPFLAGS) $< > $(SRCDIR)/$(@F)

$(OBJDIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(SRCDIR)/$(<F) -o $@

#######################################################################
#                          include rules                              #
#######################################################################

-include Makefile.inc.kernel

# vim:ft=make:ts=4:sw=4:nu:et:ai:
