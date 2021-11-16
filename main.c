#include <stdio.h>
#include <stdlib.h>

struct pixel {
	float r;
	float g;
	float b;
};

int
main(int argc, char *argv[])
{
	int width;
	int height;
	struct pixel *framebuffer;
	int i, j;
	char *filename = "output.ppm";

	width = 1024;
	height = 768;

	framebuffer = (struct pixel *)malloc(sizeof(struct pixel) * width * height);

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			framebuffer[j + i * width].r = i/(float)height;
			framebuffer[j + i * width].g = j/(float)width;
			framebuffer[j + i * width].b = 0.0;
		}
	}

	FILE *file;
	file = fopen(filename, "wb");

	if (file == NULL) {
		fprintf(stderr, "Cannot open %s for writting. Exiting.\n", filename);
		return -1;
	}

	/* write magic numbers for portable pixmap format */
	unsigned char r,g,b;
	fprintf(file, "P6\n%d %d\n255\n", width, height);
	for (i = 0; i < width * height; i++) {
		r = 255 * framebuffer[i].r;
		g = 255 * framebuffer[i].g;
		b = 255 * framebuffer[i].b;
		fwrite(&r, sizeof(char), 1, file);
		fwrite(&g, sizeof(char), 1, file);
		fwrite(&b, sizeof(char), 1, file);
	}
	fclose(file);
	free(framebuffer);
	return 0;
}
