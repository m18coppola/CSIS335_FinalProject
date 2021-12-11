#/bin/sh
make
cd bin
cd frames
rm ./*
cd ./..
./raytracer_serial
cd frames
yes | ffmpeg -framerate 24 -i %d_frame.ppm -crf 15 ../../output.mp4
