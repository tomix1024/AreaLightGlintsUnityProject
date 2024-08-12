PNG_FILES=$(wildcard RepresentativeVideo_*.png)

all : RepresentativeVideo.mp4 RepresentativeVideo_Directional.mp4

%.mp4 : %_0000.png
	ffmpeg -y -r 60 -i $(basename $@)_%04d.png -c:v libx264 -crf 20 -pix_fmt yuv420p $@
