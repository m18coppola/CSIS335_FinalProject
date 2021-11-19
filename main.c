/* 
 * Serial Ray Tracer
 *
 * Authors: Michael Coppola, Nicklaus Shelby, Andrew Towse
 * Version: Fall 2022
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cglm/cglm.h"

/* Data Types */
typedef struct {
	vec3 center;
	float radius;
} Sphere;

typedef union {
	struct {float r, g, b; };
	vec3 v;
} color_t;

/* Function Declarations */
int sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0);
void cast_ray(vec3 ray_origin, vec3 ray_direction, Sphere sphere, color_t *c);

/* Function Implementations */

/* 
 * Determines if a a interesects a sphere, and where.
 * returns 0 if there is no intersection
 * returns 1 if there is an intersection
 *
 * If there is an intersection, t0 will be populated with the distance
 * between the ray's origin and the intersection.
 *
 * See: http://www.lighthouse3d.com/tutorials/maths/ray-sphere-intersection/
 */
int
sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0)
{
	float dir_dot_v;
	float distance_from_center;
	float dist;

	/* let v be a vector from ray origin to sphere center */
	vec3 v;
	glm_vec3_sub(sphere.center, origin, v);

	/* calculate the length of the projection of the sphere center */
	dir_dot_v = glm_vec3_dot(dir, v);
	
	if (dir_dot_v > 0.0) { /* if the sphere center is in front of the ray */
		/* let pc be a vector that is the projection of the sphere center onto the ray */
		vec3 pc;
		glm_vec3_scale(dir, dir_dot_v, pc);

		/* calculate how far the center of the sphere is from the projection */
		distance_from_center = fabs(glm_vec3_distance(sphere.center, pc));

		if (distance_from_center > sphere.radius) { /* if the distance is too far from the ray, it does not intersect */
			return 0;
		} else {
			/* use pythagorean's theorem to find the leg between the radius (hypotenuse) and our distance from center (other leg) */
			dist = sqrt(sphere.radius * sphere.radius - distance_from_center * distance_from_center);

			/* the ray intersects at the projection minus the calculated leg */
			*t0 = dir_dot_v - dist;
			return 1;
			
		}
	}
	return 0;
}

void
cast_ray(vec3 ray_origin, vec3 ray_direction, Sphere sphere, color_t *color)
{
	float sphere_distance;

	sphere_distance = FLT_MAX;
	if (!sphere_ray_intersect(sphere, ray_origin, ray_direction, &sphere_distance)) {
		color->r = 0.2;
		color->g = 0.7;
		color->b = 0.8;
	} else {
		color->r = 0.4;
		color->g = 0.4;
		color->b = 0.3;
	}
}

int
main(int argc, char *argv[])
{
	int width;
	int height;
	color_t *framebuffer;
	int i, j;
	char *filename = "output.ppm";
	float fov;
	

	width = 1024;
	height = 768;
	fov = M_PI/2.0;
	Sphere ball = {.center = {-3, 0, -16}, .radius = 2};

	/* allocate an array of pixels for the image */
	framebuffer = (color_t *)malloc(sizeof(color_t) * width * height);

	/* color each pixel using the ray cast results */
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {

			#if 0
			/* Source: https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays */

			/* 
 			 * Imagine our 3d scene. To make an image of it would
 			 * be like standing in the scene and holding up a hollow picture frame.
 			 * In our program, the image file is similar to the picture frame.
 			 * To fill in the picture frame, we must calculate the path the light
 			 * travels through the picture frame and into our eye. In ray tracing, we
 			 * minimize workload by following the path backwards, so in a way we cast
 			 * a ray from our eye, through the frame and out onto the objects in the
 			 * scene. Once we see that the ray has intersected with an object, we can
 			 * conclude that the color of the light ray is the color of the object it
 			 * landed upon. We can cast a ray through each pixel of the frame to
 			 * determine what color each pixel should be.
 			 */

			/* Eye should always be 1.0 away from the picture frame */
			/* Camera should be pointed in the -z direction (right-hand rule) */
			
			/* each pixel's x/y value is between 0 and width/height (raster coords)*/
			/* we need to convert this to be between 0 and 1 (normalized device coords)*/
			/* we add 0.5 (half of a pixel) to center the ray onto the cell */

			float x_ndc = (i + 0.5)/(float)width;
			float y_ndc = (j + 0.5)/(float)height;

			/* now we map this from 0to1 to -1to1 */
			float x_frame = 2 * x_ndc - 1;
			float y_frame = 1 - 2 * y_ndc;

			/* We've made the assumption that the screen is square when we normalized the coords */
			/* We will correct this by multiplying our x by our aspect ratio */
			float aspect_ratio = (float)width/(float)height;
			x_frame *= aspect_ratio;

			/* now we will account for the field of view */
			/* We'll achieve this by growing or shrinking the frame height to the desired fov */
			/* we find by how much with trig */
			/* See: https://www.scratchapixel.com/images/upload/ray-tracing-camera/camprofile.png? */
			/* the frame and the eye are fixed at 1 distance unit apart */
			/* the ray assumes the fov is 90 degrees , so we can multiply by tan(fov/2) to scale proportionally */

			x_frame *= tan(fov/2.);

			/* and done */
			float x = x_frame;
			float y = y_frame;

			#else

			/* same thing but all in two lines */
			float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
			float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
			#endif

			vec3 dir;
			dir[0] = x;
			dir[1] = y;
			dir[2] = -1;
			glm_vec3_normalize(dir);
			vec3 origin = {0,0,0};
			cast_ray(origin, dir, ball, &(framebuffer[j + i * width]));
		}
	}

	FILE *file;
	file = fopen(filename, "wb"); /* w for write flag, b for binary mode */

	if (file == NULL) {
		fprintf(stderr, "Cannot open %s for writting. Exiting.\n", filename);
		return -1;
	}

	/* write magic numbers for portable pixmap format */
	/* P6 as magic number means full-color binary-type image file */
	/* followed by width then height */
	/* See: https://en.wikipedia.org/wiki/Netpbm */
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


