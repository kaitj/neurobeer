#!/bin/bash

# Git
DIR=/opt/git
mkdir -p $DIR/dmtcp

git clone https://github.com/dmtcp/dmtcp $DIR/dmtcp
cd $DIR/dmtcp

# Install
./configure
make 
make install 
