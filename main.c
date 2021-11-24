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
	Color color;
	float specular_reflectance;
	float diffuse_reflectance;
	float mirror_reflectance;
	float shininess;
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
const Material DARKBLUE = {.color.v = {0.2, 0.7, 0.8}};
const Material DARKRED =  {.color.v = {0.3, 0.1, 0.1} , .specular_reflectance = 0.1, .diffuse_reflectance = 0.9, .shininess = 10.0, .mirror_reflectance = 0.0};
const Material IVORY =  {.color.v = {0.4, 0.4, 0.3} , .specular_reflectance = 0.3, .diffuse_reflectance = 0.6, .shininess = 50.0, .mirror_reflectance = 0.1};
const Material MIRROR =  {.color.v = {1.0, 1.0, 1.0} , .specular_reflectance = 10.0, .diffuse_reflectance = 0.0, .shininess = 1425.0, .mirror_reflectance = 0.8};

/* Function Declarations */
int sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0);
Color cast_ray(vec3 ray_origin, vec3 ray_direction, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth);
int first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, int sphere_count, vec3 *hit, vec3 *N, Material *mat);

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
	/* let L be a vector from the ray's origin and the center of the sphere*/
	vec3 L;
	glm_vec3_sub(sphere.center, origin, L);
	
	/* let tc be the distance between the origin and */
	/* the sphere center's projection onto the ray */
	float tc = glm_vec3_dot(L, dir);

	/* let d2 be the distance between the sphere center and its projection, squared */
	float d2 = glm_vec3_dot(L, L) - tc * tc;

	float radius_squared = sphere.radius * sphere.radius;

	if (d2 > radius_squared) return 0; //sphere is too far from ray

	/* let t1c be the distance between the sphere's projection onto the ray and the point of intersection */
	float t1c = sqrt(radius_squared - d2);

	/* there could be two intersections, one from the ray entering the sphere, and one exiting the sphere */
	/* t0 is first guess, t1 is second guess */
	*t0 = tc - t1c;
	float t1 = tc + t1c;
	if (*t0 < 0) *t0 = t1; // the first guess was behind the ray origin, use the other
	if (*t0 < 0) return 0; //the second guess was also behind the ray origin, no intersection
	return 1;
}

/* 
 * Reflects the first vector over the axis the second vector represents.
 * The result is stored in the third vector
 */
void
reflect(vec3 light_dir, vec3 normal, vec3 reflection_out)
{
	glm_vec3_scale(normal, 2.0 * glm_vec3_dot(light_dir, normal), reflection_out);
	glm_vec3_sub(light_dir, reflection_out, reflection_out);
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
cast_ray(vec3 ray_origin, vec3 ray_direction, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth)
{
	vec3 intersection_pos, surface_normal;
	Material material;
	Color reflect_color;

	/* return bg color if no intersection */
	if (depth == 0 || !first_intersect_of(ray_origin, ray_direction, spheres, sphere_count, &intersection_pos, &surface_normal, &material)) {
		return DARKBLUE.color;
	}

	/* calc reflection */
	vec3 reflect_dir;
	reflect(ray_direction, surface_normal, reflect_dir);

	vec3 reflection_origin;
	vec3 perturbation;
	glm_vec3_scale(surface_normal, 1e-3, perturbation);
	if (glm_vec3_dot(reflect_dir, surface_normal) < 0.0) {
		glm_vec3_sub(intersection_pos, perturbation, reflection_origin);
	} else {
		glm_vec3_add(intersection_pos, perturbation, reflection_origin);
	}
	reflect_color = cast_ray(reflection_origin, reflect_dir, spheres, sphere_count, lights, light_count, depth - 1);

	/* calculate how much light is on point */
	float diffuse_intensity = 0;
	float specular_intensity = 0;
	for (int i = 0; i < light_count; i++) {
		/* calc light direction */
		vec3 light_dir;
		glm_vec3_sub(lights[i].pos, intersection_pos, light_dir);
		glm_vec3_normalize(light_dir);

		


		
		/* calc shadow */
		float light_distance = glm_vec3_distance(lights[i].pos, intersection_pos);

		vec3 shadow_origin;
		/* we slightly move the shadow origin from the surface of the object */
		/* to prevent it from intersecting with the sphere is lies on */
		if (glm_vec3_dot(light_dir, surface_normal) < 0.0) {
			glm_vec3_sub(intersection_pos, perturbation, shadow_origin);
		} else {
			glm_vec3_add(intersection_pos, perturbation, shadow_origin);
		}

		vec3 shadow_intersection, shadow_normal;
		Material temp_mat;
		/* we're shadowed if we intersect with a sphere before we hit the light */
		int shadowed = (first_intersect_of(shadow_origin, light_dir, spheres, sphere_count, &shadow_intersection, &shadow_normal, &temp_mat) && glm_vec3_distance(shadow_intersection, shadow_origin) < light_distance);

		if (!shadowed) {
			/* diffuse component */
			diffuse_intensity += lights[i].intensity * glm_max(0.0, glm_vec3_dot(light_dir, surface_normal));
			/* specular component */
			vec3 light_reflection;
			reflect(light_dir, surface_normal, light_reflection);
			specular_intensity += pow(glm_max(0.0, glm_vec3_dot(light_reflection, ray_direction)), material.shininess) * lights[i].intensity;
		}
	}

	/* create new color by multiplying by diffuse component, and adding the shine */
	Color output_color;
	glm_vec3_scale(material.color.v, diffuse_intensity * material.diffuse_reflectance, output_color.v);
	vec3 shine = {1.0, 1.0, 1.0};
	glm_vec3_scale(shine, specular_intensity * material.specular_reflectance, shine);
	glm_vec3_add(output_color.v, shine, output_color.v);
	glm_vec3_scale(reflect_color.v, material.mirror_reflectance, reflect_color.v);
	glm_vec3_add(output_color.v, reflect_color.v, output_color.v);
	return output_color;
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
first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, int sphere_count, vec3 *hit, vec3 *surface_normal, Material *mat)
{
	float sphere_distance = FLT_MAX;
	for (int i = 0; i < sphere_count; i++) {
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
	fov = 57;
	fov *= M_PI/180;

	Sphere spheres[] = {
		{.center = {-3.0, 0, -16}, .radius = 2, .mat = IVORY},
		{.center = {-1.0, -1.5, -12}, .radius = 2, .mat = MIRROR},
		{.center = {1.5, -0.5, -18}, .radius = 3, .mat = DARKRED},
		{.center = {7, 5, -18}, .radius = 4, .mat = MIRROR}
	};

	Light lights[] = {
		{.pos = {-20, 20, 20}, .intensity = 1.5},
		{.pos = {30, 50, -25}, .intensity = 1.8},
		{.pos = {30, 20, 30}, .intensity = 1.7}
	};

	/* allocate an array of pixels for the image */
	framebuffer = (Color *)malloc(sizeof(Color) * width * height);

	/* color each pixel using the ray cast results */
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {

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
			framebuffer[i + j * width] = cast_ray(origin, dir, spheres, 4, lights, 3, 4);
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
		Color color = framebuffer[i];

		#if 1
		/* normalize our colors */
		/* it'll dim the other channels to adjust for overly bright blow-outs in other channels */
		float max = glm_max(color.v[0], glm_max(color.v[1], color.v[2]));
		if (max > 1.0) glm_vec3_scale(color.v, 1.0 / max, color.v);
		
		#else
		/* this will just cap the brightness instead */
		/* only channels of the color that are blowing-out are dimmed */
		for (int ind = 0; ind < 3; ind++) {
			if (color.v[ind] > 1.0) {
				color.v[ind] = 1.0;
			}
		}
		#endif

		/* our rgb channels range from 0.0 to 1.0 */
		/* ppm expects 0 to 255, so we must convert */
		r = 255 * color.r;
		g = 255 * color.g;
		b = 255 * color.b;

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


