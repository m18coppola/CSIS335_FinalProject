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

	/* allocate an array of pixels for the image */
	framebuffer = (struct pixel *)malloc(sizeof(struct pixel) * width * height);

	/* set arbitrary colors for each pixel*/
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			framebuffer[j + i * width].r = i/(float)height;
			framebuffer[j + i * width].g = j/(float)width;
			framebuffer[j + i * width].b = 0.0;
		}
	}

	FILE *file;
	file = fopen(filename, "wb"); /* w for write flag, b for binary mode */

	if (file == NULL) {
		fprintf(stderr, "Cannot open %s for writting. Exiting.\n", filename);
		return -1;
	}

	/* write magic numbers for portable pixmap format */
	/* P6 means means full-color binary-type image file */
	/* followed by width then height */
	fprintf(file, "P6\n%d %d\n255\n", width, height);

	/* dump each pixel into the file */
	unsigned char r,g,b;
	for (i = 0; i < width * height; i++) {
		/* our rgb channels range from 0.0 to 1.0 */
		/* ppm expects 0 to 255, so we must convert */
		r = 255 * framebuffer[i].r;
		g = 255 * framebuffer[i].g;
		b = 255 * framebuffer[i].b;

		/* write each channel of the pixel */
		/* fwrite(adress_of_data, size_of_data_type, number_of_elements, file) */
		fwrite(&r, sizeof(char), 1, file);
		fwrite(&g, sizeof(char), 1, file);
		fwrite(&b, sizeof(char), 1, file);
	}
	fclose(file);

	/* clean up */
	free(framebuffer);

	return 0;
}
