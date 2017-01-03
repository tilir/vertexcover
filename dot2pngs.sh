#!/bin/sh

mkdir -p pngs
mkdir -p dots
mkdir -p mpss

mv *.mps mpss

for i in *.dot 
do
  dot -Tpng $i > pngs/$i.png
  mv $i dots
done

