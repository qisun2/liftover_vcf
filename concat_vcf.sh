#!/bin/bash

for i in {1..10}
do
cd input_dir$i
for j in {1..10}
do

if [  $j !=  $i ]
then 
#  File name format: tmpchr3.txt
pwd
echo adding ../input_dir${j}/tmpchr${i}.txt to tmpchr${i}.txt
cat ../input_dir${j}/tmpchr${i}.txt >> tmpchr${i}.txt
fi

done
cd ..
done
