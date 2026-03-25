#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif

static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);

double sqr(double x) { return x * x; };

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const {
		return sqrt(norm2());
	}
	void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double& operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
	Ray(const Vector& origin, const Vector& unit_direction) : O(origin), u(unit_direction) {};
	Vector O, u;
};

class Object {
public:
	Object(const Vector& albedo, bool mirror = false, bool transparent = false) : albedo(albedo), mirror(mirror), transparent(transparent) {};

	virtual bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const = 0;

	Vector albedo;
	bool mirror, transparent;
};

class Sphere : public Object {
public:
	Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent), C(center), R(radius) {};

	// returns true iif there is an intersection between the ray and the sphere
	// if there is an intersection, also computes the point of intersection P, 
	// t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
	// and the unit normal N
	bool intersect(const Ray& ray, Vector& P, double &t, Vector& N) const {
		 // TODO (lab 1) : compute the intersection (just true/false at the begining of lab 1, then P, t and N as well)
		Vector L = ray.O - C;
		double b = 2 * dot(ray.u, L);
		double c = dot(L, L) - R * R;
		//self - note: check the definition on the slides once again, don't rely on your goldfish memory alone.
		double delta = b * b - 4 * c;
		if (delta < 0) return false;
		double t1 = (-b - sqrt(delta)) / 2.0; //check if t1 or t2 is positive or not too (note from the lecture.)
		//smallest is t1,  from the thing but we want the smallest positive one.
		//distance to the camera initially remember too.
		double t2 = (-b + sqrt(delta)) / 2.0; //check with the professor on this, not too sure on this.
		if (t2 < 0) return false;
		//We want the smallest possible t
		t = (t1 >= 0) ? t1 : t2;
		P = ray.O + t * ray.u; // P(t) = O + t*u
		N = P - C;
		N.normalize();
		return true;
	}

	double R;
	Vector C;
};


// I will provide you with an obj mesh loader (labs 3 and 4)
class TriangleMesh : public Object {
public:
	TriangleMesh(const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent) {};

	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
		// TODO (labs 3 and 4)
		return false;
	}
};


class Scene {
public:
	Scene() {};
	void addObject(const Object* obj) {
		objects.push_back(obj);
	}

	// returns true iif there is an intersection between the ray and any object in the scene
    // if there is an intersection, also computes the point of the *nearest* intersection P, 
    // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
    // and the unit normal N. 
	// Also returns the index of the object within the std::vector objects in object_id
	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N, int &object_id) const  {

		// TODO (lab 1): iterate through the objects and check the intersections with all of them, 
		// and keep the closest intersection, i.e., the one if smallest positive value of t
		t = 1e9; // Initializing to a very large value (infinity) I guess.
		bool has_intersection = false;
		//Idea is to loop through all the intersections to find the nearest thing.
		for (int i = 0; i < objects.size(); i++) {
			Vector P_cur, N_cur; //self note- debug this might break during compilation.
			double t_cur;
			if (objects[i]->intersect(ray, P_cur, t_cur, N_cur)) {
				//basicly we found the intersection, is it closer we do that here.
				if (t_cur < t) {
					t = t_cur;
					P = P_cur;
					N = N_cur;
					object_id = i;
					has_intersection = true;
				}
			}
		}
		return has_intersection;

	}


	// return the radiance (color) along ray
	Vector getColor(const Ray& ray, int recursion_depth) {

		if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);

		// TODO (lab 1) : if intersect with ray, use the returned information to compute the color ; otherwise black 
		// in lab 1, the color only includes direct lighting with shadows

		Vector P, N;
		double t;
		int object_id;
		if (intersect(ray, P, t, N, object_id)) {

			if (objects[object_id]->mirror) {
				Vector R = ray.u - 2.0 * dot(ray.u, N) * N;
				// Vector R = ray.u - 2.0 * dot(ray.u, N-1) * N;
				Ray R_ray(P + 1e-4 * N, R);
				// return getColor in the reflected direction, with recursion_depth+1 (recursively)
				return getColor(R_ray, recursion_depth + 1);
			} // else

			if (objects[object_id]->transparent) { // optional

				// return getColor in the refraction direction, with recursion_depth+1 (recursively)
			} // else

			// test if there is a shadow by sending a new ray
			// if there is no shadow, compute the formula with dot products etc.
			Vector L = light_position - P;
			double d = L.norm();
			L.normalize();

			//Imp: Shadow ray, offset slightly to avoid self-intersection
			Ray shadow_ray(P + 1e-4 * N, L);
			Vector P_s, N_s;
			double t_s;
			int object_id_s;

			//you also wanna change the value and stuff to make sure things are considered.
			// If we hit an object and it's closer than the light (t_s < d), we are in shadow!
			if (intersect(shadow_ray, P_s, t_s, N_s, object_id_s) && t_s < d) {
				return Vector(0, 0, 0); // Pitch black shadow I assume?
			}
			double intensity = light_intensity / (4.0 * M_PI * d * d);
			double diffuse_factor = std::max(0.0, dot(N, L)); // Based on what we were taught in the lecture that we don't consider the negative stuff.

			// Note: Sir, I panicked a bit and couldn't recall this properly, according to me - Color = Intensity * Cos(theta) * (Albedo / Pi) // but as per notes eq: L = (I / 4*pi*d^2) * (rho / pi) * max(0, <N, w_i>), so implemented what seemed more compliant to me.
			Vector diffuse_color = objects[object_id]->albedo * (intensity * diffuse_factor / M_PI);

			return diffuse_color;

			// TODO (lab 2) : add indirect lighting component with a recursive call
		}



		return Vector(0, 0, 0);
	}

	std::vector<const Object*> objects;

	Vector camera_center, light_position;
	double fov, gamma, light_intensity;
	int max_light_bounce;
};


int main() {
	int W = 512;
	int H = 512;

	for (int i = 0; i<32; i++) {
		engine[i].seed(i);
	}

	Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8));
	Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
	Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
	Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
	Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));

	Scene scene;
	scene.camera_center = Vector(0, 0, 55); //Change in lecture.
	scene.light_position = Vector(-10,20,40);
	scene.light_intensity = 3E7;
	scene.fov = 60 * M_PI / 180.;
	scene.gamma = 1.0;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.max_light_bounce = 5;

	scene.addObject(&center_sphere);

	/*
	scene.addObject(&wall_left);
	scene.addObject(&wall_right);
	scene.addObject(&wall_front);
	scene.addObject(&wall_behind);
	scene.addObject(&ceiling);
	scene.addObject(&floor);
	*/

	std::vector<unsigned char> image(W * H * 3, 0);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector color;

			// TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)			
			double x = j - W / 2.0 + 0.5; //Center horizontal pixel mapping.
			double y = H / 2.0 - i - 0.5; //Center the vertical mapping (basically invert the y-axis like explained int he class.)
			double z = -W / (2.0 * tan(scene.fov / 2.0)); //Map focal lenght via W and FOV.
			//New Ray and new origin remember.
			
			// Vector ray_direction(0., 0., -1);
			Vector ray_direction(x, y, z);
			//Adding the ray direction normalization here too.
			ray_direction.normalize(); //we must ensure that the lenght is exactly 1.


			Ray ray(scene.camera_center, ray_direction);

			// TODO (lab 2) : add Monte Carlo / averaging of random ray contributions here
			// TODO (lab 2) : add antialiasing by altering the ray_direction here
			// TODO (lab 2) : add depth of field effect by altering the ray origin (and direction) here

			color  = scene.getColor(ray, 0);

			image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(color[0] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(color[1] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(color[2] / 255., 1. / scene.gamma)));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}

//Notes for self during lecture -
//Don't forget the object id and all too.
//Put white pixel if there is an intersection for example.
//don't forget to return the object_id.
//if image is totally white it might be possible your camera is inside the thing.