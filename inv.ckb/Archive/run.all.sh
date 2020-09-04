#!/bin/bash
if [ $# -lt 6 ]
	then
		echo './run.all.sh xc yc zc xwd ywd zwd'
		exit 1
fi
xc=$1
yc=$2
zc=$3
xwd=$4
ywd=$5
zwd=$6

# edit the following script before excuting run.all.csh

./inv.ckb.make.sh ${xc} ${yc} ${zc} ${xwd} ${ywd} ${zwd}
./inv.ckb.joint.sh ${xwd} ${ywd} ${zwd}
./run.ckb.syn.sh ${xwd} ${ywd} ${zwd}
./run.solver.1th.sh ${xwd} ${ywd} ${zwd}
