#!/bin/bash
if [ $# -lt 6 ]
	then
		echo './inv.ckb.make.sh xc yc zc xwd ywd zwd'
		exit 1
fi
xc=$1
yc=$2
zc=$3
xwd=$4
ywd=$5
zwd=$6

bin_make_ckb='../../../codes/inv.ckb/inv_make_ckb'

# dimensions of the inversion blocks
mx=160; my=106; mz=25;
echo $mx $my $mz
# checkerboard pattern
# 1 = linear square
# 2 = linear cylinder
# 3 = linear spher
# 4 = cos^2 cylinder
# 5 = cos^2 spher
itype=1

# index of first ckb block
nx1=1; ny1=1; nz1=1
# index of last ckb block
nx2=$mx; ny2=$my; nz2=$mz

# maximum perturbation
dvmax=0.10

# Checkerboard with 5x5 inversion cells in lat and lon and 1-cell edge
mxc=${xc}; wid1=${xwd}; edgex1=0; edgex2=0; dvmin1=$dvmax
myc=${yc}; wid2=${ywd}; edgey1=0; edgey2=0; dvmin2=$dvmax
mzc=${zc}; wid3=${zwd}; edgez1=0; edgez2=0; dvmin3=$dvmax

${bin_make_ckb} << EOF
$itype
$mx $my $mz
$mxc $myc $mzc
$dvmax
$wid1 $edgex1 $edgex2 $dvmin1
$wid2 $edgey1 $edgey2 $dvmin2
$wid3 $edgez1 $edgez2 $dvmin3
$nx1 $ny1 $nz1
$nx2 $ny2 $nz2
EOF

cp block_ckb.dat 'block_ckb_'${xwd}'x'${ywd}'x'${zwd}'.dat'
