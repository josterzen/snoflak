#!/bin/sh

[ -d build ] || mkdir build
[ -d dist ] || mkdir dist
[ -d dist/frames ] || mkdir dist/frames
[ -d dist/figures ] || mkdir dist/figures
cd build

cmake ..
cmake --build .

mv ./snowflake ../dist/snowflake
cd ..
