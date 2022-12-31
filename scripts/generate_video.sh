#!/bin/sh

cd ../dist/figures

ffmpeg -framerate 60 -pattern_type glob -i '*.png' -vf scale=1920:-2,setsar=1:1 -c:v libx265 video.mp4
