#!/bin/bash
# this script will compile the file specified in the command line argument
# if no file specified it will compile all .cpp files in the src/ directory

# --- Variables
CXX=g++
SRC_PATH=src
SRC_EXT=cpp
EXE_EXT=exe
COMPILE_FLAGS="-std=c++17 -Wall -Wextra -g -O0"
INCLUDES="-I include/"

# --- Command Line Arguments
if [ "$#" -eq 1 ]; then # if there is one argument
    source_file=$1

    # remove directories or endings
    src_file=${source_file%.$SRC_EXT}
    src_file=${src_file%.$EXE_EXT}
    src_file=${src_file#*$SRC_PATH/}
    echo $src_file

    # compile
    $CXX $COMPILE_FLAGS $INCLUDES include/funcs.cpp $SRC_PATH/$src_file.$SRC_EXT -o $src_file.$EXE_EXT

else # otherwise to all files
    for source_file in $SRC_PATH/*.$SRC_EXT;
    do
        # extract filename
        src_file=${source_file%.$SRC_EXT}
        src_file=${src_file#*$SRC_PATH/}
        echo $src_file
        
        # compile
        $CXX $COMPILE_FLAGS $INCLUDES include/funcs.cpp $source_file -o $src_file.$EXE_EXT
        echo "-----------"
    done
fi

exit 0