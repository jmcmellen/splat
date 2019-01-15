#!/bin/sh

for i in `ls *.sdf 2>/dev/null`
do
    file=`echo $i | tr ':' '_'`
    echo "mv $i $file"
    mv $i $file
done
for i in `ls *.sdf.bz2 2>/dev/null`
do
    file=`echo $i | tr ':' '_'`
    echo "mv $i $file"
    mv $i $file
done
