#/bin/sh
gmake
cd bin
cd frames
rm ./*
cd ./..
#./raytracer_serial
mpiexec -n 2 ./raytracer_threaded
cd frames
yes | ffmpeg -framerate 24 -i %d_frame.ppm -crf 15 ../../output.mp4
