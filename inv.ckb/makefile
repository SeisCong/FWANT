#******************************************************************************#
#*                         makefile of LSQR solver                            *#
#******************************************************************************#

FC=ifort
all:
	${FC} -o inv_ckb_integral -O0 -g -static-intel inv_ckb_integral.f90
	${FC} -o inv_ckb_synthetic -O0 -g -static-intel inv_ckb_synthetic.f90
	${FC} -o inv_make_ckb -O0 -g -static-intel inv_make_ckb.f90

# vim:ft=make:ts=4:sw=4:nu:et:ai:
