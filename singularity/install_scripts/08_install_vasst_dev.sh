#!/bin/bash

# Install vasst-dev repo & MCR v92 (2017a)
DIR=/opt/git/vasst-dev
HASH=9cfced6a089a857adcf1e3a2e709c4917f16c592

# MCR
TMP_DIR=/tmp/mcr
MCR_DIR=$DIR/mcr/v92
mkdir -p $MCR_DIR $TMP_DIR
MCR_DIR=`realpath $MCR_DIR`

curl -L --retry 5 http://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip >> $TMP_DIR/install.zip
pushd $TMP_DIR
unzip install.zip
./install -mode silent -agreeToLicense yes -destinationFolder $MCR_DIR
popd

rm -rf $TMP_DIR

# VASST-DEV
git clone https://github.com/akhanf/vasst-dev $DIR/
pushd $DIR
git checkout $HASH
popd

# if octave exists:
# if [ -e /etc/octave.conf ]
# then
#     echo addpath\(genpath\(\'${PIPELINE_TOOL_DIR}/matlab\'\)\)\; >> /etc/octave.conf
# fi
