#!/bin/sh
# This wrapper script is intended to support independent execution.
# 
# This script uses the following environment variables set by the submit MATLAB code:
# MDCE_MATLAB_EXE     - the MATLAB executable to use
# MDCE_MATLAB_ARGS    - the MATLAB args to use
#

# Copyright 2010-2017 The MathWorks, Inc.

# Only saw this with the Hydra version, but to avoid this, will load a later version of zlib
# 
#    /share/pkg/mpich2/1.4.1p1/bin/mpiexec "/share/pkg/matlab/R2017a/bin/worker" -parallel
#    /share/pkg/mpich2/1.4.1p1/bin/mpiexec: /lib64/libz.so.1: version `ZLIB_1.2.3.3' not found
#                          (required by /share/pkg/matlab/R2017a/bin/glnxa64/libxml2.so.2) 
#
module load zlib/1.2.3.9

export TZ="America/New_York"
echo "Executing: ${MDCE_MATLAB_EXE} ${MDCE_MATLAB_ARGS}"
eval "${MDCE_MATLAB_EXE}" ${MDCE_MATLAB_ARGS}
EXIT_CODE=${?}
echo "Exiting with code: ${EXIT_CODE}"
exit ${EXIT_CODE}
