#!/bin/bash

DIR=/opt/git/
mkdir -p $DIR/neurobeer

# Git
git clone https://github.com/kaitj/neurobeer $DIR/neurobeer
cd $DIR/neurobeer

# Install requirements
pip3 install -r requirements.txt

# Install neurobeer
python3 setup.py install
