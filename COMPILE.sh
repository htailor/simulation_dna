#!/bin/sh
./clear_all.sh
EXE_FILE=Nucleation
cp $EXE_FILE "$EXE_FILE.bak" 
cd ./src
clear
echo "Compiling new Nucelation file..."
echo "==============================================================================="
make -f Makefile
echo "==============================================================================="
echo "COMPLETE."
rm *.o
mv $EXE_FILE ..
cd ..
echo
