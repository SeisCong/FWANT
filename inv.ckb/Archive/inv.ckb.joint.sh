#!/bin/bash
# you may assign different Vp and Vs ckeckerboard
if [ $# -lt 3 ]
	then
		echo './inv.ckb.joint.sh xwd ywd zwd'
		exit 1
fi
xwd=$1
ywd=$2
zwd=$3

ckb_Vp='block_ckb.dat'
ckb_Vs='block_ckb.dat'

# and change Vp and Vs magnitudes with fct
fct_Vp=1.0
fct_Vs=1.0

ckb_out='block_ckb_all_'${xwd}'x'${ywd}'x'${zwd}'.dat'

nVp=`cat $ckb_Vp|wc -l`
nVs=`cat $ckb_Vs|wc -l`

#cat /dev/null > $ckb_out
echo "$nVp + $nVs" | bc > $ckb_out

awk '{print $1*"'$fct_Vp'"}' $ckb_Vp >> $ckb_out
awk '{print $1*"'$fct_Vs'"}' $ckb_Vs >> $ckb_out
