#!/bin/bash

if [[ "$#" -ne 2 ]]; then
    echo "Illegal number of parameters"
fi

echo "the directory is"
echo $1

echo "comment is "
echo $2

cd $1

git commit -m $2 -a




