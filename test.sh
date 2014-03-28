#!/bin/sh

echo === VERIFY TOOL BUILDING ===
make -C tools
if [ ! $? -eq 0 ]; then
	echo Building failed with error $?.
	exit
fi
echo == VERIFY TOOL BUILT ===

echo === COMPILING START ===
make
if [ ! $? -eq 0 ]; then
	echo Compiling failed with error $?.
	exit
fi
echo === COMPILING END ===

echo === EXECUTION START ===
./3d2rewq data/input_test.in input_test.out input_test.log
echo === EXECUTION END ===

echo === DIFF START ===
tools/verify input_test.out data/output.correct
echo === DIFF END ===

