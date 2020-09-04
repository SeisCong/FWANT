#!/bin/bash

outdir="result.2th"
fnm_smooth='../block.49x45x26.1x1x1.2x2x1.smooth2th.dat'

# input list file name
list_Gd="inv_Gd_list"

# damping list
# ite_01
#damp_list="0.5 1 2 4 8 16 32 64 128 256"
# ite_02
damp_list="0.5 1 2 4 8 16 32 64 128 256"

# smoothness list
# ite_01
#smot_list="1 4 16 64 256"
# ite_02
smot_list="0 1 4 16 64"


#weigst_list="20"
weigst_list="0"
#weigev_list="50"
weigev_list="0"
#weiglo_list="80"
weiglo_list="0"

# station terms per data
nst=1
# event terms per data
nev=1
# location terms per data
nlo=3

# num_cmp should be 2 when both Vp and Vs are inverted
num_cmp=2

# see ./sim.kernel/run.kernel.assm.sh for ker_thres and figures from ./inv.*/inv_plot_G.m
#ker_thres="1e-2"
ker_thres="1e-3"

# absolute phase delay less than
data_thres=20

#SOLVER='../../code/inv.LSQR/solver_damp_smooth'
SOLVER='/net/fs01/data/yang/n.cascadia/codes/inv.LSQR/solver_damp_smo_loc'

if [ ! -d $outdir ]; then
   mkdir -p $outdir
fi

# ----------------------------------------------------------------------
function inv_compile {
  #FC=pgf90
  FC=ifort
  SRC='/net/fs01/data/yang/n.cascadia/codes/inv.LSQR'
  echo "$FC -O3 -o $SOLVER $SRC/mod_LSQR.f90 $SRC/solver_damp_smooth.f90"
  $FC -O3 -o $SOLVER $SRC/mod_LSQR.f90 $SRC/solver_damp_smooth.f90
}

# ----------------------------------------------------------------------
function inv_loop {
  for weig_lo in ${weiglo_list}; do
  for weig_st in ${weigst_list}; do
  for weig_ev in ${weigev_list}; do
  for smot in ${smot_list}; do
  for damp in ${damp_list}; do
  
  echo "loc st ev smooth damp=" $weig_lo $weig_st $weig_ev $smot $damp

$SOLVER << EOF
1
$list_Gd
$fnm_smooth
$num_cmp
$nst $weig_st
$nev $weig_ev
$nlo $weig_lo
$damp $smot
$ker_thres $data_thres
EOF

  mv try.xyz $outdir/try.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.sum $outdir/sum.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.err $outdir/err.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.sta $outdir/sta.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.evt $outdir/evt.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  mv try.loc $outdir/loc.damp${damp}.smot${smot}.st${weig_st}.ev${weig_ev}.lo${weig_lo}.dat
  
  done
  done
  done
  done
  done
}

# ----------------------------------------------------------------------
#echo "compiling the program ..."
#time inv_compile;

echo ""
echo "solving with different daming ..."
time inv_loop;
