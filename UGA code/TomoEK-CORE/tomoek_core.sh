#!/bin/bash
HOME="/home/labm/.core"
homefolder='/home/labm/.core/TomoEK'
includefolder='include'
srcfolder='src'
nodeid=`echo "$(pwd)" | grep -Eo "[[:digit:]]+" | tail -n1`

clear

ifconfig eth0 broadcast 172.16.255.255 
/usr/lib/quagga/zebra -u root -g root
/usr/lib/quagga/batmand eth0
ip rule add from all lookup 66
subfolder=$nodeid

sleep 10

clear
mkdir $subfolder
#mkdir $includefolder
cp $homefolder/$subfolder/list.dat $subfolder 
$homefolder/data/ >> datapath.dat 
#cp -a $homefolder/include/. $includefolder
#cp -a $homefolder/src/. $srcfolder  
#g++ $homefolder/hello.cpp
#sudo cp /home/labm/.core/TomoEK/half1.dat ~/ >> move.dat
#./a.out
#chmod 777 TomoEK
#./TomoEK >> output.dat
$HOME/TomoEK/TomoEK $homefolder $subfolder/list.dat $subfolder>> out.dat
ls >> final.dat


