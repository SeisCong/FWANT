#!/bin/bash
if [ $# -lt 3 ]
	then
		echo './run.ckb.syn.sh xwd ywd zwd'
		exit 1
fi
xwd=$1
ywd=$2
zwd=$3

bin_ckb_syn='../../../codes/inv.ckb/inv_ckb_synthetic'

# input list file name

list_Gd_rcd="inv_Gd_list_observed"

fnm_ckb='block_ckb_all_'${xwd}'x'${ywd}'x'${zwd}'.dat'
list_Gd='inv_Gd_list_'${xwd}'x'${ywd}'x'${zwd}'.ori'

# num_cmp should be 2 when both Vp and Vs are inverted
num_cmp=2
# see ./run.kernel.assm.sh in sim.kernel
ker_thres="5e-4"

# maximum absolute synthetic delay time
data_thres=10

# ----------------------------------------------------------------------
${bin_ckb_syn} << EOF
$num_cmp
$fnm_ckb
$list_Gd_rcd
$list_Gd
$ker_thres $data_thres
EOF
