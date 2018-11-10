#!/bin/bash

rm -rf $1
mkdir $1

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o $1/ldstat -s -O2 -std=c++11 -static

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o $1/ldstat.exe -s -O2 -std=c++11 -static

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o $1/ldstat.exe -s -O2 -std=c++11 -static

elif [ $1 == "macos" ]; then

    LIB=/usr/local/lib

    LIBGCC=${LIB}/gcc/x86_64-apple-darwin17.5.0/8.1.0/libgcc.a

    g++ *.cpp -o $1/ldstat -O2 -std=c++11 ${LIBGCC}

fi
