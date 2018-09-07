#!/bin/bash

# Install vasst-dev repo
DIR=/opt/git/vasst-dev

# VASST-DEV
git clone https://github.com/kaitj/vasst-dev $DIR/

echo addpath\(genpath\(\'$DIR/tools/matlab\'\)\)\; >> /etc/octave.conf
