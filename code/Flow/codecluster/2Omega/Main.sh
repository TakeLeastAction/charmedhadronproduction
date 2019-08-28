#!/bin/bash

curlist=$(pwd)
echo $curlist
for dirlist in $(ls $curlist)
do
echo $dirlist
if test -d $dirlist
then

cp makejob ./${dirlist}
cp makejob-condor ./${dirlist}
cp run.sh ./${dirlist}
cd ${dirlist}
./makejob
  curlist2=$(pwd) 
  for dir1 in $(ls $curlist2)
  do
    if test -d $dir1
     then
       cp input.txt ./${dir1}
       cp circle.txt ./${dir1}
       cp jobNum.txt ./${dir1}
       cp PT.txt ./${dir1}

       
    fi
  done
./makejob-condor
cd ..

fi
done


