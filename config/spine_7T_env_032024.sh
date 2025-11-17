#!/bin/bash

toolbox_home=/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox
anaconda_dir=/export02/data/landelle/anaconda/

# matlab
LD_PREFIX="/export01/local/matlab23b/sys/os/glnxa64:/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/libraries"
export  LD_LIBRARY_PATH=/export01/local/matlab23b/bin/glnxa64/



export TMPDIR='/export02/data/tmp'


export PATH="/cerebro/cerebro1/dataset/spine_7T/derivatives/toolboxes/spinalcordtoolbox/bin:$PATH"
export SCT_DIR=/cerebro/cerebro1/dataset/spine_7T/derivatives/toolboxes/spinalcordtoolbox/
export MPLBACKEND=Agg

# python setup

source /export02/data/landelle/anaconda/etc/profile.d/conda.sh
conda activate CL_spine_7T_env_py9
echo "++ Python version adjusted : `which python`"

# ANTS setup
export PATH=${PATH}:/usr/lib/ants


# FSL Setup
FSLDIR=/cerebro/cerebro1/dataset/bmpd/derivatives/thibault_test/code/toolbox/fsl
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
. ${FSLDIR}/etc/fslconf/fsl.sh


# AFNI configuration
PREFIX=$toolbox_home/afni
if [ -e $PREFIX ] ; then
   if [ -z "$PATH" ] ; then
      export PATH="$PREFIX"
   else
      export PATH="$PREFIX:$PATH"
   fi
   echo "++ AFNI added to PATH: $PREFIX"
fi




