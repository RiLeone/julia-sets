#!/usr/bin/env bash
ffmpeg -r 10 -i plts/frames/%05d.png -vcodec libx264 -pix_fmt yuv420p plts/vids/zooming_experience.mp4
ffmpeg -r 10 -i plts/frames/%05d.png plts/gifs/julia_animation.gif
#rm -rf plts/frames/*.png
