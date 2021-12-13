/* 
 * Serial Software Ray Tracer
 *
 * Authors: Michael Coppola
 * Version: Fall 2022
 */
#include <cglm/cglm.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

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

/* Function Declarations */
Color cast_ray(vec3 orig, vec3 dir, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth);
int first_intersect_of(vec3 origin, vec3 dir, Sphere *spheres, int sphere_count, vec3 *hit, vec3 *N, Material *mat);
void refract(vec3 incomming, vec3 normal, float refractive_index, vec3 out);
void render(Color *framebuffer, int width, int height, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth, double fov);
int sphere_ray_intersect(Sphere sphere, vec3 origin, vec3 dir, float *t0);
void reflect(vec3 light_dir, vec3 normal, vec3 reflection_out);
int write_frame(Color *framebuffer, int width, int height, long unsigned int frameID);

/* Material Definitions */
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
 * See: http://www.lighthouse3d.com/tutorials/maths/ray-sphere-intersection/
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
 * param origin ray origin
 * param dir direction of ray
 * param spheres array of spheres to be considered for intersection
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

void
render(Color *framebuffer, int width, int height, Sphere *spheres, int sphere_count, Light *lights, int light_count, int depth, double fov)
{
	int i, j;
	float x, y;
	for (j = 0; j < height; j++)
		for (i = 0; i < width; i++) {
			x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
			y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
			vec3 dir;
			dir[0] = x;
			dir[1] = y;
			dir[2] = -1;
			glm_vec3_normalize(dir);
			vec3 origin = {0,0,0};
			framebuffer[i + j * width] = cast_ray(origin, dir, spheres, sphere_count, lights, light_count, depth);
		}
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

	int initReturn, rank, numProcs;
	initReturn = MPI_Init(NULL, NULL);
	if(initReturn != MPI_SUCCESS) {
	  fprintf(stderr, "MPI was unable to initialize. Exiting.\n");
	  return 1;
	}

	
	/* get info from MPI context */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	
	/* rotation speed */
	degs_per_sec = 90.0;
	rads_per_sec = degs_per_sec * M_PI/180.0;

	duration = 5.0; //seconds
	framerate = 24; //frames per second

	frame_count = framerate * duration;
	time_step = duration / (float)frame_count;
	

	width = 720;
	height = 480;
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
	int ready_thread;
	framebuffer = (Color *)malloc(sizeof(Color) * width * height);
	
	  // IF RANK 0 wait for proc to request for work and hand it out
	  if (rank == 0){
	    // Loop to send proc work
	    
	    for(int frame = 0; frame < frame_count; frame++) {
	      // send to waiting process frame number to calculate
	      printf("Receiving which thread to send work to\n");
	      MPI_Recv(&ready_thread, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      printf("Sending frame to %d\n", ready_thread);
	      MPI_Send(&frame, 1, MPI_INT, ready_thread, 1, MPI_COMM_WORLD);  
	    }
	    
	    // send -1 to each process as frame to do
	    int frame = -1;
	    for(int procs = 0; procs < numProcs-1; procs++){
	      MPI_Recv(&ready_thread, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      MPI_Send(&frame, 1, MPI_INT, ready_thread, 1, MPI_COMM_WORLD);
	      printf("[%d] Finished\n", ready_thread);
	    }
	  }
	  else{
	    int recv_frame = 0;
	    while(recv_frame != -1){
	      // Sending to rank 0 we are ready
	      printf("[%d] Sending Ready to rank 0\n", rank);
	      MPI_Send(&rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	      
	      // Request Work from rank 0
      	      MPI_Recv(&recv_frame, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      printf("[%d] Received Frame %d\n",rank, recv_frame);

	      // IF frame to do  == -1 end loop
	      if(recv_frame == -1){
		break;
	      }
	      // Adjust every sphere in the frame
	      	for (int si = 0; si < sphere_count; si++) {
			spheres[si].center[2] += 20.0;
			glm_vec3_rotate(spheres[si].center, rads_per_sec * (time_step * recv_frame), GLM_YUP);
			spheres[si].center[2] -= 20.0;
			
		}
		
		render(framebuffer, width, height, spheres, sphere_count, lights, light_count, reflective_depth, fov);
		if (write_frame(framebuffer, width, height, recv_frame) == -1) {
		  return -1;
		}

		// Undo rotations
	       	for (int si = 0; si < sphere_count; si++) {
		        spheres[si].center[2] += 20.0;
			glm_vec3_rotate(spheres[si].center, -rads_per_sec * (time_step * recv_frame), GLM_YUP);
			spheres[si].center[2] -= 20.0;
		}
	
	    }
      	}

	  // MPI_Finalize
	  MPI_Finalize();
       

	/* clean up */
	free(framebuffer);

	return 0;
}


