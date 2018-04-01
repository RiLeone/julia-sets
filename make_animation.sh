#!/usr/bin/env bash
ffmpeg -r 10 -i plts/frames/%05d.png plts/gifs/julia_animation.gif
rm -rf plts/frames/*.png
