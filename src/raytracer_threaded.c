/* 
 * Serial Software Ray Tracer
 *
 * Authors: Michael Coppola, Andrew Towse, Nick Shelby
 * Version: Fall 2022
 */
#include <cglm/cglm.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

/* Data Types */

typedef union {
	struct {float r, g, b; };
	float arr[3];
	vec3 v;
} Color;

typedef struct {
	vec3 pos;
	float intensity;
} Light;

typedef struct {
	Color color;

	float diffuse_reflectance;
	float specular_reflectance;
	float mirror_reflectance;
	float refraction_reflectance;

	float shininess;
	float refractive_index;
} Material;

typedef struct {
	vec3 center;
	float radius;
	Material mat;
} Sphere;

typedef struct {
	Color *framebuffer;
	int width;
	int height;
	Sphere *spheres;
	int sphere_count;
	Light *lights;
	int light_count;
	int depth;
	float fov;
} SharedData;

typedef struct {
	pthread_t thread_id;
	SharedData *td;
	int start_row;
	int end_row;
} Thread;

/* Function Declarations */
Color cast_ray(vec3 orig, vec3 dir, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth);
int first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, int sphere_count, vec3 *hit, vec3 *N, Material *mat);
void reflect(vec3 light_dir, vec3 normal, vec3 reflection_out);
void refract(vec3 incomming, vec3 normal, float refractive_index, vec3 out);
void render(Color *framebuffer, int width, int height, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth, float fov, int thread_count);
int sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0);
void * thread_cast(void *args);
int write_frame(Color *framebuffer, int width, int height, long unsigned int frameID);

/* Material Definitions */
/* 
 * DARKBLUE !INCOMPLETE MATERIAL - COLOR ONLY!
 * DARKRED
 * WHITE
 * YELLOW
 * IVORY
 * MIRROR
 * GLASS
 */
const Material DARKBLUE = {.color.v = {0.2, 0.7, 0.8}};
const Material DARKRED =  {
	.color.v = {0.3, 0.1, 0.1} ,
	.diffuse_reflectance = 0.9, 
	.specular_reflectance = 0.1,
	.mirror_reflectance = 0.0,
	.refraction_reflectance = 0.0,
	.shininess = 10.0, 
	.refractive_index = 1.0,
};
const Material WHITE =  {
	.color.v = {0.3, 0.3, 0.3} ,
	.diffuse_reflectance = 0.9, 
	.specular_reflectance = 0.1,
	.mirror_reflectance = 0.0,
	.refraction_reflectance = 0.0,
	.shininess = 10.0, 
	.refractive_index = 1.0,
};
const Material YELLOW =  {
	.color.v = {0.3, 0.2, 0.1} ,
	.diffuse_reflectance = 0.9, 
	.specular_reflectance = 0.1,
	.mirror_reflectance = 0.0,
	.refraction_reflectance = 0.0,
	.shininess = 10.0, 
	.refractive_index = 1.0,
};
const Material IVORY =  {
	.color.v = {0.4, 0.4, 0.3} , 
	.diffuse_reflectance = 0.6, 
	.specular_reflectance = 0.3, 
	.mirror_reflectance = 0.1,
	.refraction_reflectance = 0.0,
	.shininess = 50.0, 
	.refractive_index = 1.0,
};
const Material MIRROR =  {
	.color.v = {1.0, 1.0, 1.0} , 
	.diffuse_reflectance = 0.0, 
	.specular_reflectance = 10.0, 
	.mirror_reflectance = 0.8,
	.refraction_reflectance = 0.0,
	.shininess = 1425.0, 
	.refractive_index = 1.0,
};
const Material GLASS =  {
	.color.v = {0.6, 0.7, 0.8} , 
	.diffuse_reflectance = 0.0, 
	.specular_reflectance = 0.5, 
	.mirror_reflectance = 0.1,
	.refraction_reflectance = 0.8,
	.shininess = 125.0, 
	.refractive_index = 1.5,
};

/* Function Implementations */

/* 
 * casts a ray of virtual light into the world backwards and
 * returns the color of it after it interacts with the scene
 *
 * param orig: origin location of the ray
 * param dir: direction of the ray
 * param spheres: list of spheres in scene
 * param sphere_count: sphere list length
 * param lights: list of light sources in scene
 * param light_count: light list length
 * param depth: number of allowed reflections a ray can take off a surface
 *
 * return: color of light ray
 */
Color
cast_ray(vec3 orig, vec3 dir, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth)
{
	vec3 hit, surface_normal;
	Material material;
	Color reflect_color;
	Color refract_color;
	vec3 perturbation;
	vec3 reflect_dir;
	vec3 refract_dir;
	vec3 reflection_orig;
	vec3 refraction_orig;
	int i;
	float diffuse_intensity = 0;
	float specular_intensity = 0;
	float light_dist;
	vec3 light_dir;
	vec3 shadow_origin;
	vec3 shadow_intersection;
	vec3 shadow_normal;
	Material temp_mat;
	int shadowed;
	vec3 light_reflection;
	vec3 shine = {1.0, 1.0, 1.0};
	vec3 refraction_component;
	Color output_color;

	/* background */
	if (depth == 0 || !first_intersect_of(orig, dir, spheres, sphere_count, &hit, &surface_normal, &material))
		return DARKBLUE.color;

	/* get reflection & refraction */
	glm_vec3_scale(surface_normal, 1e-3, perturbation);
	reflect(dir, surface_normal, reflect_dir);
	refract(dir, surface_normal, material.refractive_index, refract_dir);
	if (glm_vec3_dot(reflect_dir, surface_normal) < 0.0)
		glm_vec3_sub(hit, perturbation, reflection_orig);
	else
		glm_vec3_add(hit, perturbation, reflection_orig);
	if (glm_vec3_dot(refract_dir, surface_normal) < 0.0)
		glm_vec3_sub(hit, perturbation, refraction_orig);
	else
		glm_vec3_add(hit, perturbation, refraction_orig);
	reflect_color = cast_ray(reflection_orig, reflect_dir, spheres, sphere_count, lights, light_count, depth - 1);
	refract_color = cast_ray(refraction_orig, refract_dir, spheres, sphere_count, lights, light_count, depth - 1);

	for (i = 0; i < light_count; i++) {
		/* determine if there is a shadow */
		glm_vec3_sub(lights[i].pos, hit, light_dir);
		glm_vec3_normalize(light_dir);
		light_dist = glm_vec3_distance(lights[i].pos, hit);
		if (glm_vec3_dot(light_dir, surface_normal) < 0.0)
			glm_vec3_sub(hit, perturbation, shadow_origin);
		else
			glm_vec3_add(hit, perturbation, shadow_origin);
		shadowed = (first_intersect_of(shadow_origin, light_dir, spheres, sphere_count, &shadow_intersection, &shadow_normal, &temp_mat) && glm_vec3_distance(shadow_intersection, shadow_origin) < light_dist);

		/* illuminate surface */
		if (!shadowed) {
			diffuse_intensity += lights[i].intensity * glm_max(0.0, glm_vec3_dot(light_dir, surface_normal));
			reflect(light_dir, surface_normal, light_reflection);
			specular_intensity += pow(glm_max(0.0, glm_vec3_dot(light_reflection, dir)), material.shininess) * lights[i].intensity;
		}
	}

	/* blend diffuse, specular, reflective and refractive color channels */
	glm_vec3_scale(material.color.v, diffuse_intensity * material.diffuse_reflectance, output_color.v);
	glm_vec3_scale(shine, specular_intensity * material.specular_reflectance, shine);
	glm_vec3_add(output_color.v, shine, output_color.v);
	glm_vec3_scale(reflect_color.v, material.mirror_reflectance, reflect_color.v);
	glm_vec3_add(output_color.v, reflect_color.v, output_color.v);
	glm_vec3_scale(refract_color.v, material.refraction_reflectance, refraction_component);
	glm_vec3_add(refraction_component, output_color.v, output_color.v);
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
 * param origin: ray origin
 * param dir: direction of ray
 * param spheres: array of spheres to be considered for intersection
 * output hit: coordinate of intersection
 * output surface_normal: unit vector perpendicular to surface of intersection
 * output material: material of surface intersected
 *
 * return: 1 if intersection, 0 otherwise
 */

int
first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, int sphere_count, vec3 *hit, vec3 *surface_normal, Material *mat)
{
	float distance = FLT_MAX;
	float plane_dist = FLT_MAX;
	float current_distance;
	float d;
	vec3 pt;
	int i;

	for (i = 0; i < sphere_count; i++)
		if (sphere_ray_intersect(spheres[i], origin, dir, &current_distance) && current_distance < distance) {
			distance = current_distance;
			glm_vec3_scale(dir, current_distance, *hit);
			glm_vec3_add(origin, *hit, *hit);
			glm_vec3_sub(*hit, spheres[i].center, *surface_normal);
			glm_vec3_normalize(*surface_normal);
			*mat = spheres[i].mat;
		}

	/* CHECKERBOARD MATH */
	//TODO: factor this out somehow
	if (fabs(dir[1]) > 1e-3) {
		d = -(origin[1]+4)/dir[1];
		glm_vec3_scale(dir, d, pt);
		glm_vec3_add(origin, pt, pt);
		if (d > 0 && fabs(pt[0]) < 10 && pt[2] < -10 && pt[2] > -30 && d < distance) {
			plane_dist = d;
			glm_vec3_dup(pt, *hit);
			glm_vec3_zero(*surface_normal);
			(*surface_normal)[1] = 1.0;
			if (((int)(0.5 * (*hit)[0] + 1000) + (int)(0.5 * (*hit)[2])) & 1)
				*mat = WHITE;
			else
				*mat = YELLOW;
		}
	}

	return glm_min(distance, plane_dist) < 1000;
}
/* 
 * renders the input spheres, lights and parameters to the supplied framebuffer.
 *
 * output framebuffer: framebuffer to write the render output to
 * param width: width of frame buffer
 * param height: height of framebuffer
 * param spheres: list of spheres in scene
 * param sphere_count: sphere list length
 * param lights: list of light sources in scene
 * param light_count: light list length
 * param depth: upper limit to number of reflections a ray can make
 * param fov: the fov of the perspective of the frame (in radians)
 */
void
render(Color *framebuffer, int width, int height, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth, float fov, int thread_count)
{
	int current_row;
	int i;
	int rc;
	int row_count;

	SharedData sd = {
		.framebuffer = framebuffer,
		.width = width,
		.height = height,
		.spheres = spheres,
		.sphere_count = sphere_count,
		.lights = lights,
		.light_count = light_count,
		.depth = depth,
		.fov = fov,
	};

	Thread *threads = (Thread *)malloc(sizeof(Thread) * thread_count);
	current_row = 0;
	row_count = height/thread_count;
	for (i = 0; i < thread_count; i++) {
		threads[i].td = &sd;
		threads[i].start_row = current_row;
		current_row += row_count;
		if (i < height%thread_count) current_row++;
		threads[i].end_row = current_row;
		rc = pthread_create(&(threads[i].thread_id), NULL, thread_cast, threads + i);
		if (rc != 0) {
			fprintf(stderr, "Could not create thread %d.\n", i);
			exit(1);
		}
	}
	for (i = 0; i < thread_count; i++)
		pthread_join(threads[i].thread_id, NULL);
	free(threads);
}

/* 
 * Reflects the first vector over the the same axis as the normal vector's.
 * The reflection is stored in reflection_out
 */
void
reflect(vec3 light_dir, vec3 normal, vec3 reflection_out)
{
	glm_vec3_scale(normal, 2.0 * glm_vec3_dot(light_dir, normal), reflection_out);
	glm_vec3_sub(light_dir, reflection_out, reflection_out);
}

/* 
 * creates a new vector that is a refraction of the incoming vector.
 * The incomming angle is calculated with the incomming vector and the surface normal
 * the outgoing vector is calculated with the refractive index, and stored in the out vector
 *
 * This function uses textbook Snell's law.
 * See: https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
 */
void
refract(vec3 incomming, vec3 normal, float refractive_index, vec3 out)
{
	float cosi;
	float etai = 1.0;
	float etat, eta;
	vec3 n;
	float k;
	float temp;

	cosi = -1.0 * glm_max(-1.0, glm_min(1.0, glm_vec3_dot(incomming, normal)));
	etat = refractive_index;
	glm_vec3_dup(normal, n);
	if (cosi < 0) {
		cosi = -cosi;
		temp = etai;
		etai = etat;
		etat = temp;
		glm_vec3_negate(n);
	}
	eta = etai / etat;
	k = 1 - eta * eta * (1 - cosi * cosi);
	if (k < 0)
		glm_vec3_zero(out);
	else {
		glm_vec3_scale(incomming, eta, out);
		glm_vec3_scale(n, eta * cosi - sqrt(k), n);
		glm_vec3_add(out, n, out);
	}
}

/* 
 * Determines if a ray interesects a sphere, and where.
 * returns 0 if there is no intersection
 * returns 1 if there is an intersection
 *
 * If there is an intersection, t0 will be populated with the distance
 * between the ray's origin and the intersection.
 *
 * See: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
 */
int
sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0)
{
	vec3 L;
	float tc;
	float d2;
	float r2;
	float t1c;
	float t1;

	glm_vec3_sub(sphere.center, origin, L);
	tc = glm_vec3_dot(L, dir);
	d2 = glm_vec3_dot(L, L) - tc * tc;
	r2 = sphere.radius * sphere.radius;
	if (d2 > r2) return 0;
	t1c = sqrt(r2 - d2);
	*t0 = tc - t1c;
	t1 = tc + t1c;
	if (*t0 < 0) *t0 = t1;
	if (*t0 < 0) return 0;
	return 1;
}

/* 
 * thread execution function for render.
 *
 * param args: a fully populated Thread pointer
 */
void *
thread_cast(void *args)
{
	Thread this = *(Thread*)args;
	SharedData ta = *(this.td);

	int i, j;
	float x, y;
	for (j = this.start_row; j < this.end_row; j++)
		for (i = 0; i < ta.width; i++) {
			x =  (2*(i + 0.5)/(float)ta.width  - 1)*tan(ta.fov/2.)*ta.width/(float)ta.height;
			y = -(2*(j + 0.5)/(float)ta.height - 1)*tan(ta.fov/2.);
			vec3 dir;
			dir[0] = x;
			dir[1] = y;
			dir[2] = -1;
			glm_vec3_normalize(dir);
			vec3 origin = {0,0,0};
			ta.framebuffer[i + j * ta.width] = cast_ray(origin, dir, ta.spheres, ta.sphere_count, ta.lights, ta.light_count, ta.depth);
		}
	return NULL;
}

/* 
 * write the framebuffer to a ppm file in the ./frame/ sub-directory.
 *
 * param framebuffer: framebuffer to print to file
 * param width: width of framebuffer
 * param height: height of framebuffer
 * param frameID: UID for output filename
 */
int
write_frame(Color *framebuffer, int width, int height, long unsigned int frameID)
{
	int i;
	char filename[32];
	FILE *file;
	unsigned char r,g,b;
	float max;
	Color color;

	sprintf(filename, "frames/%lu_frame.ppm", frameID);
	file = fopen(filename, "wb");
	if (file == NULL) {
		fprintf(stderr, "Cannot open %s for writting. Exiting.\n", filename);
		return -1;
	}
	/* write magic numbers for portable pixmap format */
	/* P6 as magic number means full-color binary-type image file */
	/* followed by width then height */
	/* See: https://en.wikipedia.org/wiki/Netpbm */
	fprintf(file, "P6\n%d %d\n255\n", width, height);
	for (i = 0; i < width * height; i++) {
		color = framebuffer[i];
		max = glm_max(color.v[0], glm_max(color.v[1], color.v[2]));
		if (max > 1.0) glm_vec3_scale(color.v, 1.0 / max, color.v);
		r = 255 * color.r;
		g = 255 * color.g;
		b = 255 * color.b;
		fwrite(&r, sizeof(char), 1, file);
		fwrite(&g, sizeof(char), 1, file);
		fwrite(&b, sizeof(char), 1, file);
	}
	fclose(file);
	return 0;
}

int
main(int argc, char *argv[])
{
	int width;
	int height;
	Color *framebuffer;
	int reflective_depth;
	float fov;
	float duration;
	int framerate;
	long unsigned int frame_count;
	float time_step;
	float degs_per_sec;
	float rads_per_sec;
	int thread_count;
	char ffmpeg_command[256];
	int i;

	/* rotation speed */
	degs_per_sec = 90.0;
	rads_per_sec = degs_per_sec * M_PI/180.0;

	duration = 5.0; //seconds
	framerate = 24; //frames per second

	frame_count = framerate * duration;
	time_step = duration / (float)frame_count;

	width = 400;
	height = 300;
	reflective_depth = 4;
	fov = 57;
	fov *= M_PI/180;

	Sphere spheres[] = {
		{.center = {-3.0, 0, -16}, .radius = 2, .mat = IVORY},
		{.center = {-1.0, -1.5, -12}, .radius = 2, .mat = GLASS},
		{.center = {1.5, -0.5, -18}, .radius = 3, .mat = DARKRED},
		{.center = {7, 5, -18}, .radius = 4, .mat = MIRROR}
	};
	int sphere_count = 4;

	Light lights[] = {
		{.pos = {-20, 20, 20}, .intensity = 1.5},
		{.pos = {30, 50, -25}, .intensity = 1.8},
		{.pos = {30, 20, 30}, .intensity = 1.7}
	};
	int light_count = 3;

	thread_count = 4;

	system("exec rm -rd ./frames/");
	system("exec mkdir ./frames/");

	framebuffer = (Color *)malloc(sizeof(Color) * width * height);
	for (int frame = 0; frame < frame_count; frame++) {
		render(framebuffer, width, height, spheres, sphere_count, lights, light_count, reflective_depth, fov, thread_count);

		for (int si = 0; si < sphere_count; si++) {
			spheres[si].center[2] += 20.0;
			glm_vec3_rotate(spheres[si].center, rads_per_sec * time_step, GLM_YUP);
			spheres[si].center[2] -= 20.0;
			
			if (write_frame(framebuffer, width, height, frame) == -1) {
				return -1;
			}
		}
		printf("\r%.0f%% of frames rendered.", (float)frame/(float)frame_count*100.0);
		fflush(stdout);

	}
	printf("\r100%% of frames rendered.");
	fflush(stdout);
	printf("\nRendering complete.\n");
	printf("elapsed time: ??\n");
	printf("Handing off to ffmpeg...\n");
	printf("[");
	/* super secret "goes to" operator ;) */
	i = 68; while (i-->0) printf("=");
	printf("]\n");
	sprintf(ffmpeg_command, "exec ffmpeg -framerate %d -i ./frames/%%d_frame.ppm -crf 15 -y output.mp4", framerate);
	printf("%s\n", ffmpeg_command);
	system(ffmpeg_command);

	/* clean up */
	free(framebuffer);

	return 0;
}
