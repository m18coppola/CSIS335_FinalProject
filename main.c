/* 
 * Serial Software Ray Tracer
 *
 * Authors: Michael Coppola
 * Version: Fall 2022
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cglm/cglm.h"

/* Data Types */

typedef union {
	struct {float r, g, b; };
	vec3 v;
} Color;

typedef struct {
	Color diffuse_color;
	// other members will be added later
} Material;

typedef struct {
	vec3 center;
	float radius;
	Material mat;
} Sphere;

typedef struct {
	vec3 pos;
	float intensity;
} Light;

/* Constants */
const Material DARKBLUE = {.diffuse_color.v = {0.2, 0.7, 0.8} };
const Material DARKRED =  {.diffuse_color.v = {0.3, 0.1, 0.1} };
const Material IVORY =  {.diffuse_color.v = {0.4, 0.4, 0.3} };

/* Function Declarations */
int sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0);
Color cast_ray(vec3 ray_origin, vec3 ray_direction, Sphere *spheres, Light *lights);
int first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, vec3 *hit, vec3 *N, Material *mat);

/* Function Implementations */

/* 
 * Determines if a ray interesects a sphere, and where.
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
		distance_from_center = glm_vec3_distance(sphere.center, pc); 

		/* if the distance between the ray and the center is less than the radius, no intersection */
		if (distance_from_center > sphere.radius) {
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

/* 
 * Takes a ray and traces it through space until it intersects with an object
 *
 * returns the material of the object it lands on or the material of the background if
 * there is no intersection
 *
 * param ray_origin the starting position of the ray
 * param ray_direction unit vector representing the direction of the ray
 * param spheres array of spheres to be considered for intersection testing
 * param lights array of lights to shine on spheres
 */
Color
cast_ray(vec3 ray_origin, vec3 ray_direction, Sphere *spheres, Light *lights)
{
	vec3 point, surface_normal;
	Material material;

	/* return bg color if no intersection */
	if (!first_intersect_of(ray_origin, ray_direction, spheres, &point, &surface_normal, &material)) {
		return DARKBLUE.diffuse_color;
	}

	/* otherwise, calculate how much light is on it */
	/* for each light, check angle between surface normal and light direction */
	/* use that angle as a parameter to adjust brightness of output pixel */
	/* the further the angle of the surface is tilted away from the light, proportionally darken */
	float diffuse_light_intensity = 0;
	for (int i = 0; i < 1; i++) { //TODO: fix hard coded range
		/* vector subtraction then normalize to get unit vector */
		vec3 light_dir;
		glm_vec3_sub(lights[i].pos, point, light_dir);
		glm_vec3_normalize(light_dir);
		/* get angle by using the dot product of light direction and surface normal*/
		float angle = glm_vec3_dot(light_dir, surface_normal);
		/* add to the brightness along with the other lights */
		diffuse_light_intensity += lights[i].intensity * ((angle > 0) ? angle : 0); 
	}
	/* create new color by multiplying the original one with the intensity scalar */
	Color illuminated_color;
	glm_vec3_scale(material.diffuse_color.v, diffuse_light_intensity, illuminated_color.v);
	return illuminated_color;
}

/*
 * get the closest intersection of the input ray and the input spheres
 * returns 1 if there is an intersection, and 0 otherwise.
 * 
 * if there is an intersection, this function populates the parameter pointers with:
 *
 * *hit the position of the intersection
 * *surface_normal a unit vector that is perpendicular with the tangential plane upon the position of intersection
 * *mat the material of the intersected object
 *
 * param origin ray origin
 * param dir direction of ray
 * param spheres array of spheres to be considered for intersection
 */

int
first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, vec3 *hit, vec3 *surface_normal, Material *mat)
{
	float sphere_distance = FLT_MAX;
	for (int i = 0; i < 4; i++) { //TODO: fix hard coded range
		float current_distance;
		/* if this intersection is closer than the last calculated one */
		if (sphere_ray_intersect(spheres[i], origin, dir, &current_distance) && current_distance < sphere_distance) {
			/* update distance */
			sphere_distance = current_distance;

			/* the product of the distance on the direction is the position */
			glm_vec3_scale(dir, current_distance, *hit);
			/* transform from world coords to camera coords */
			glm_vec3_add(origin, *hit, *hit);
			/* make surface normal with the difference between the sphere center and the intersection */
			glm_vec3_sub(*hit, spheres[i].center, *surface_normal);
			glm_vec3_normalize(*surface_normal);
			/* get material */
			*mat = spheres[i].mat;
		}
	}
	/* 1000 is the hardcoded clipping distance */
	/* should meet our needs for now */
	return sphere_distance < 1000;
}

int
main(int argc, char *argv[])
{
	int width;
	int height;
	Color *framebuffer;
	int i, j;
	char *filename = "output.ppm";
	float fov;
	

	width = 1024;
	height = 768;
	fov = 90;
	fov *= M_PI/180;

	Sphere spheres[] = {
		{.center = {7, 5, -18}, .radius = 4, .mat = IVORY},
		{.center = {1.5, -0.5, -18}, .radius = 3, .mat = DARKRED},
		{.center = {-1.0, -1.5, -12}, .radius = 2, .mat = DARKRED},
		{.center = {-3, 0, -16}, .radius = 2, .mat = IVORY}
	};

	Light lights[] = {
		{.pos = {-20, 20, 20}, .intensity = 1.5}
	};

	/* allocate an array of pixels for the image */
	framebuffer = (Color *)malloc(sizeof(Color) * width * height);

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
			framebuffer[j + i * width] = cast_ray(origin, dir, spheres, lights);
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


