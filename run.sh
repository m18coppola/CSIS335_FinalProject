#/bin/sh
make
cd bin
./raytracer_serial
cd frames
ffmpeg -framerate 24 -i %d_frame.ppm -crf 15 ../../output.mp4
