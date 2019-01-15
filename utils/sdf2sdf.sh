#!/bin/sh

for i in `ls *.sdf`
do
    file=`echo $i | tr ':' '_'`
    echo "mv $i $file"
    mv $i $file
done
for i in `ls *.sdf.bz2`
do
    file=`echo $i | tr ':' '_'`
    echo "mv $i $file"
    mv $i $file
done
