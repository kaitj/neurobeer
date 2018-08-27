#!/bin/bash

# Install vasst-dev repo
DIR=/opt/git/vasst-dev
HASH=9cfced6a089a857adcf1e3a2e709c4917f16c592

git clone https://github.com/akhanf/vasst-dev $DIR/
pushd $DIR
git checkout $HASH
popd

# if octave exists:
if [ -e /etc/octave.conf ]
then
    echo addpath\(genpath\(\'${PIPELINE_TOOL_DIR}/matlab\'\)\)\; >> /etc/octave.conf
fi
