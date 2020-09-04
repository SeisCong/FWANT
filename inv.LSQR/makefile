#******************************************************************************#
#*                         makefile of LSQR solver                            *#
#******************************************************************************#

FC=ifort #/opt/intel/composer_xe_2013_sp1.2.139/bin/intel64/ifort
FFLAG=-O0 -static -intel

#FC=pgf90
#FFLAG=-O3 -Bstatic

all:
#	${FC} -o solver_damp_smooth ${FFLAG} mod_LSQR.f90 solver_damp_smooth.f90
	${FC} -g -o solver_damp_smo_loc_yang ${FFLAG} mod_LSQR.f90 solver_damp_smo_loc_yang.f90
#	${FC} -o solver_damp_smo_loc2 ${FFLAG} mod_LSQR.f90 solver_damp_smo_loc2.f90

# vim:ft=make:ts=4:sw=4:nu:et:ai:
