#!/bin/bash

# Install dependencies
apt-get install -y openjdk-9-jdk-headless

# Install camino
DIR=/opt/camino

if [ -d $DIR ]; then
    echo Removing dir...
    rm -rf $DIR
fi

mkdir -p $DIR

# Release from: https://sourceforge.net/projects/camino
HASH=a844005b39c2fc7fd5d33c6af977e44248de42cd
VERSION=2018.08.01
URL=https://downloads.sourceforge.net/project/camino/camino-code-$HASH.zip

# In MBs -- make larger if Java virtual memory needed
# HEAPSIZE=16000

# Install camino
echo -n "Install Camino $VERSION..."
curl -L --retry 5 $URL -o $DIR/camino-code-$HASH.zip

pushd $DIR
unzip camino-code-$HASH.zip
mv camino-code-$HASH/* .
make
popd

# Update profile
# PROFILE=~/.bashrc
#
# if grep -xq "PATH=$DIR/bin:\$PATH" $PROFILE # return - if exists
# then
#     echo "PATH=$DIR/bin:" in the PATH already
# else
#     echo "export PATH=$DIR/bin:\$PATH" >> $PROFILE
#     echo "export LD_LIBRARY_PATH=$DIR/lib:/$LD_LIBRARY_PATH" >> $PROFILE
#     echo "export MANPATH=$DIR/man:\$MANPATH" >> $PROFILE
#     echo "export CAMINO_HEAP_SIZE=$HEAPSIZE" >> $PROFILE
# fi
#
# source $PROFILE
