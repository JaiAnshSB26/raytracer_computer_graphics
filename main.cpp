#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
//Imports for Lab2
#include <omp.h>
#include <algorithm>
//New Imports done for Lab2.

//Imports for LAB 3
#include <map>
#include <string>
#include <fstream>
//New imports done for Lab 3
//Needed for some Lab 3/4 methods - 
#include <list>

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
Vector operator*(const Vector& a, const Vector& b){
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]); // Element-wise multiplication needed for indirect lighting
}

//Things needed for Lab 2 defined on the very top.
//Firstly, the random cosine-weighted vector around a normal (Box-Muller adaptation.)
Vector random_cos(const Vector& N, int thread_id){
	double r1 = uniform(engine[thread_id]);
	double r2 = uniform(engine[thread_id]);

	double x = cos(2.0 * M_PI * r1) * sqrt(1.0 - r2); 
	double y = sin(2.0 * M_PI * r1) * sqrt(1.0 - r2);
	double z = sqrt(r2);

	//We now need to generate the orthogonal tangent frame for the normal mapping.
	Vector T1;
	if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) {
		T1 = Vector(0, -N[2], N[1]);
	} else if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])) {
		T1 = Vector(-N[2], 0, N[0]);
	} else {
		T1 = Vector(-N[1], N[0], 0);
	}
	T1.normalize();
	Vector T2 = cross(N, T1);

	//Now as per the lectures, we need to transform the local Box-Muller point into the world space relative to the Normal.
	return x * T1 + y * T2 + z * N;
}

//Now we need to ensure the proper uniform random distribution across threads using Box-Muller.
void boxMuller(double stdev, double& x, double& y, int thread_id){
	double r1 = uniform(engine[thread_id]);
	double r2 = uniform(engine[thread_id]);
	x = sqrt(-2.0 * log(r1)) * cos(2.0 * M_PI * r2) * stdev;
	y = sqrt(-2.0 * log(r1)) * sin(2.0 * M_PI * r2) * stdev;
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
		//distance to the camera initially remember too. (lecture note for 2nd td.)
		double t2 = (-b + sqrt(delta)) / 2.0; //check with the professor on this, not too sure on this.
		if (t2 < 1e-5) return false;
		// //We want the smallest possible t
		// if (t2 < 0) return false; // Both intersections are behind the camera

		// // We want the smallest positive t
		// t = (t1 >= 0) ? t1 : t2;
		t = (t1 > 1e-5) ? t1 : t2; //Should I use 1e-5 or 0?
		P = ray.O + t * ray.u; // P(t) = O + t*u
		N = P - C;
		N.normalize();
		return true;
	}

	double R;
	Vector C;
};


//Implementation of new classes needed for Lab 3 - I begin with this while understanding the Traingles class.
class BoundingBox {
public:
	BoundingBox() : Bmin(Vector(1e9, 1e9, 1e9)), Bmax(Vector(-1e9, -1e9, -1e9)) {}
	BoundingBox(const Vector& Bmin, const Vector& Bmax) : Bmin(Bmin), Bmax(Bmax) {}
	Vector Bmin, Bmax;

	bool intersect(const Ray& ray, double& t_inter) const {
		double tx0 = (Bmin[0] - ray.O[0]) / ray.u[0];
		double tx1 = (Bmax[0] - ray.O[0]) / ray.u[0];
		if (tx0 > tx1) std::swap(tx0, tx1);

		double ty0 = (Bmin[1] - ray.O[1]) / ray.u[1];
		double ty1 = (Bmax[1] - ray.O[1]) / ray.u[1];
		if (ty0 > ty1) std::swap(ty0, ty1);

		double tz0 = (Bmin[2] - ray.O[2]) / ray.u[2];
		double tz1 = (Bmax[2] - ray.O[2]) / ray.u[2];
		if (tz0 > tz1) std::swap(tz0, tz1);

		double t0 = std::max(tx0, std::max(ty0, tz0));
		double t1 = std::min(tx1, std::min(ty1, tz1));

		if (t0 <= t1 && t1 > 0) {
            t_inter = t0 > 0 ? t0 : t1;
            return true;
        }
        return false;
	}
};

//Even though we will be implementing the BVHNode class in Lab 4, I am defining it here, now itself for better versioning and readability.
class BVHNode {
public:
	BoundingBox bbox;
	int starting_triangle;
	int ending_triangle;
	BVHNode* child_left;
	BVHNode* child_right;

	BVHNode() : child_left(nullptr), child_right(nullptr) {}
};

// Class only used in labs 3 and 4 
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1) {
		vtx[0] = vtxi; vtx[1] = vtxj; vtx[2] = vtxk;
		uv[0] = uvi; uv[1] = uvj; uv[2] = uvk;
		n[0] = ni; n[1] = nj; n[2] = nk;
		this->group = group;
	};
	int vtx[3]; // indices within the vertex coordinates array
	int uv[3];  // indices within the uv coordinates array
	int n[3];   // indices within the normals array
	int group;  // face group
};

// Class only used in labs 3 and 4 
class TriangleMesh : public Object {
public:
	TriangleMesh(const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent) {
        root = new BVHNode();
    };

	// first scale and then translate the current object
	void scale_translate(double s, const Vector& t) {
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * s + t;
		}
	}

	// read an .obj file
	void readOBJ(const char* obj) {
		std::ifstream f(obj);
		if (!f) return;

		std::map<std::string, int> mtls;
		int curGroup = -1, maxGroup = -1;

		// OBJ indices are 1-based and can be negative (relative), this normalizes them
		auto resolveIdx = [](int i, int size) {
			return i < 0 ? size + i : i - 1;
		};

		auto setFaceVerts = [&](TriangleIndices& t, int i0, int i1, int i2) {
			t.vtx[0] = resolveIdx(i0, vertices.size());
			t.vtx[1] = resolveIdx(i1, vertices.size());
			t.vtx[2] = resolveIdx(i2, vertices.size());
		};
		auto setFaceUVs = [&](TriangleIndices& t, int j0, int j1, int j2) {
			t.uv[0] = resolveIdx(j0, uvs.size());
			t.uv[1] = resolveIdx(j1, uvs.size());
			t.uv[2] = resolveIdx(j2, uvs.size());
		};
		auto setFaceNormals = [&](TriangleIndices& t, int k0, int k1, int k2) {
			t.n[0] = resolveIdx(k0, normals.size());
			t.n[1] = resolveIdx(k1, normals.size());
			t.n[2] = resolveIdx(k2, normals.size());
		};

		std::string line;
		while (std::getline(f, line)) {
			// Trim trailing whitespace
			line.erase(line.find_last_not_of(" \r\t\n") + 1);
			if (line.empty()) continue;

			const char* s = line.c_str();

			if (line.rfind("usemtl ", 0) == 0) {
				std::string matname = line.substr(7);
				auto result = mtls.emplace(matname, maxGroup + 1);
				if (result.second) {
					curGroup = ++maxGroup;
				} else {
					curGroup = result.first->second;
				}
			} else if (line.rfind("vn ", 0) == 0) {
				Vector v;
				sscanf(s, "vn %lf %lf %lf", &v[0], &v[1], &v[2]);
				normals.push_back(v);
			} else if (line.rfind("vt ", 0) == 0) {
				Vector v;
				sscanf(s, "vt %lf %lf", &v[0], &v[1]);
				uvs.push_back(v);
			} else if (line.rfind("v ", 0) == 0) {
				Vector pos, col;
				if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &col[0], &col[1], &col[2]) == 6) {
					for (int i = 0; i < 3; i++) col[i] = std::min(1.0, std::max(0.0, col[i]));
					vertexcolors.push_back(col);
				} else {
					sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
				}
				vertices.push_back(pos);
			}
			else if (line[0] == 'f') {
				int i[4], j[4], k[4], offset, nn;
				const char* cur = s + 1;
				TriangleIndices t;
				t.group = curGroup;

				// Try each face format: v/vt/vn, v/vt, v//vn, v
				if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0], &j[0], &k[0], &i[1], &j[1], &k[1], &i[2], &j[2], &k[2], &offset)) == 9) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0], &j[0], &i[1], &j[1], &i[2], &j[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceUVs(t, j[0], j[1], j[2]);
				} else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0], &k[0], &i[1], &k[1], &i[2], &k[2], &offset)) == 6) {
					setFaceVerts(t, i[0], i[1], i[2]); 
					setFaceNormals(t, k[0], k[1], k[2]);
				} else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2], &offset)) == 3) {
					setFaceVerts(t, i[0], i[1], i[2]);
				}
				else continue;

				indices.push_back(t);
				cur += offset;

				// Fan triangulation for polygon faces (4+ vertices)
				while (*cur && *cur != '\n') {
					TriangleIndices t2;
					t2.group = curGroup;
					if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3], &offset)) == 3) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceUVs(t2, j[0], j[2], j[3]);
					} else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3], &offset)) == 2) {
						setFaceVerts(t2, i[0], i[2], i[3]); 
						setFaceNormals(t2, k[0], k[2], k[3]);
					} else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) == 1) {
						setFaceVerts(t2, i[0], i[2], i[3]);
					} else { 
						cur++; 
						continue; 
					}

					indices.push_back(t2);
					cur += offset;
					i[2] = i[3]; j[2] = j[3]; k[2] = k[3];
				}
			}
		}
	}

	void init_bvh() {   //Again not needed for now but I am defining it here for the sake of completeness.
		build_bvh(root, 0, indices.size());
	}

	bool TriangleIntersect(const Ray& ray, Vector& P, double& t, Vector& N, int triangle_id) const {
		TriangleIndices t_indices = indices[triangle_id];
		Vector A = vertices[t_indices.vtx[0]];
		Vector B = vertices[t_indices.vtx[1]];
		Vector C = vertices[t_indices.vtx[2]];
		//take the differences?
		Vector e1 = B - A;
		Vector e2 = C - A;
		Vector N_unormalized = cross(e1, e2);
		//IMP note to self - double, don't use float by mistake here again, it would blow up the code during compilation again and you won't be able to locate the error!
		double beta = dot(e2, cross(A - ray.O, ray.u)) / dot(ray.u, N_unormalized);
		double gamma = -dot(e1, cross(A - ray.O, ray.u)) / dot(ray.u, N_unormalized);
		double alpha = 1 - beta - gamma;

		//if (alpha < 0 || alpha > 1 && beta < 0 || beta > 1 && gamma < 0 || gamma > 1) return false;
		if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1) return false;

		double t_cur = dot(A - ray.O, N_unormalized) / dot(ray.u, N_unormalized);
		if (t_cur < 1e-5) return false;
		
		t = t_cur;
		P = ray.O + t * ray.u;
		
		if (normals.size() > 0) {
			Vector normalA = normals[t_indices.n[0]];
			Vector normalB = normals[t_indices.n[1]];
			Vector normalC = normals[t_indices.n[2]];
			N = alpha * normalA + beta * normalB + gamma * normalC;
			N.normalize(); //Note for self: Don't forget to normalize...
		} else {
			N = N_unormalized;
			N.normalize();
		}
		
		return true;
	}

	BoundingBox compute_bbox(int starting_triangle, int ending_triangle) const {
		BoundingBox b;
		for (int i = starting_triangle; i < ending_triangle; i++){
			TriangleIndices t_indices = indices[i];
			Vector A = vertices[t_indices.vtx[0]];
			Vector B = vertices[t_indices.vtx[1]];
			Vector C = vertices[t_indices.vtx[2]];

			//The mins...
			b.Bmin[0] = std::min(b.Bmin[0], std::min(A[0], std::min(B[0], C[0])));
			b.Bmin[1] = std::min(b.Bmin[1], std::min(A[1], std::min(B[1], C[1])));
			b.Bmin[2] = std::min(b.Bmin[2], std::min(A[2], std::min(B[2], C[2])));

			//Now the maxs...
			b.Bmax[0] = std::max(b.Bmax[0], std::max(A[0], std::max(B[0], C[0])));
			b.Bmax[1] = std::max(b.Bmax[1], std::max(A[1], std::max(B[1], C[1])));
			b.Bmax[2] = std::max(b.Bmax[2], std::max(A[2], std::max(B[2], C[2])));
		}
		return b;
	}

	Vector compute_barycenter(TriangleIndices t_indices) const {
		Vector A = vertices[t_indices.vtx[0]];
		Vector B = vertices[t_indices.vtx[1]];
		Vector C = vertices[t_indices.vtx[2]];
		return (A + B + C) / 3.0;  //simple logic, don't complicate the implementations for urself!
	}

	//Now I implement the get_longest method, which I would need in the main intersect function.

	int get_longest(Vector length) const {
		if (length[0] > length[1] && length[0] > length[2]) return 0;
		if (length[1] > length[0] && length[1] > length[2]) return 1;
		return 2; //Note for self: make sure to not confuse and with bitwise and in a rush!
	}

	//Method for lab 4 I guess, but I kept running into compilation errors, so I thought maybe including this and changing all relevant function attributes at once might be a good idea!
	void build_bvh(BVHNode* node, int starting_triangle, int ending_triangle) {
		node->bbox = compute_bbox(starting_triangle, ending_triangle);
		node->starting_triangle = starting_triangle;
		node->ending_triangle = ending_triangle;

		Vector diag = node->bbox.Bmax - node->bbox.Bmin;
		Vector middle_diag = node->bbox.Bmin + diag * 0.5;
		int longest_axis = get_longest(diag);

		int pivot_index = starting_triangle;
		for (int i = starting_triangle; i < ending_triangle; i++) { //we have to loop over and calculate it for every isntance.
			Vector barycenter = compute_barycenter(indices[i]);
			if (barycenter[longest_axis] < middle_diag[longest_axis]) {
				std::swap(indices[i], indices[pivot_index]);
				pivot_index++;
			}
		}
        //if (pivot_index >= starting_triangle && pivot_index >= ending_triangle && ending_triangle - starting_triangle < 8) return;
		if (pivot_index <= starting_triangle || pivot_index >= ending_triangle || ending_triangle - starting_triangle < 5) return;

		node->child_left = new BVHNode();
		node->child_right = new BVHNode();
		
		build_bvh(node->child_left, starting_triangle, pivot_index); //Recursive function as was instructed for the 3 rd to-do in the descriptions for Lab 3 and Lab 4!
		build_bvh(node->child_right, pivot_index, ending_triangle);
	}


	bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
		// lab 3 : for each triangle, compute the ray-triangle intersection with Moller-Trumbore algorithm - implemented above and used here.
		// lab 3 : once done, speed it up by first checking against the mesh bounding box - implemented above and used here.
		// lab 4 (I would do it in lab 3 since I might mess up later.): recursively apply the bounding-box test from a BVH datastructure - for the sake of completeness of the problem and me not messing up in the future I have implemented tits and bits of it right now as well.
		double temp_t;
		if (!root->bbox.intersect(ray, temp_t)) return false;

		std::list<BVHNode*> nodes_to_visit;
		nodes_to_visit.push_front(root);

		double best_inter_distance = 1e9;
		bool found_intersection = false;

		while (!nodes_to_visit.empty()) {
			BVHNode* curNode = nodes_to_visit.back();
			nodes_to_visit.pop_back();

			if (curNode->child_left) {
				double inter_distance;
				if (curNode->child_left->bbox.intersect(ray, inter_distance)) {
					if (inter_distance < best_inter_distance) {
						nodes_to_visit.push_back(curNode->child_left);  //I made bit of mess out of push_back method, some text on StackOverflow really helped me here.
					}
				}
				if (curNode->child_right->bbox.intersect(ray, inter_distance)) {
					if (inter_distance < best_inter_distance) {
						nodes_to_visit.push_back(curNode->child_right);
					}
				}
			} else {
				for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++) {
					Vector P_cur, N_cur;
					double t_cur;
					if (TriangleIntersect(ray, P_cur, t_cur, N_cur, i)) {
						if (t_cur < best_inter_distance) {
							best_inter_distance = t_cur;
							P = P_cur;
							N = N_cur;
							t = t_cur;
							found_intersection = true;
						}
					}
				}
			}
		}
		return found_intersection;
		//return false;
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	BVHNode* root;
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
		t = 1e9; // Initializing to a very large value (infinity), the general appraoch.
		bool has_intersection = false;
		//Idea is to loop through all the intersections to find the nearest thing (as taught in the lecture).
		for (int i = 0; i < objects.size(); i++) {
			Vector P_cur, N_cur; //self note- debug this might break during compilation, check this block again.
			double t_cur;
			if (objects[i]->intersect(ray, P_cur, t_cur, N_cur)) {
				//basicly we found the intersection, is it closer we do that check/computatioon here.
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
	Vector getColor(const Ray& ray, int recursion_depth, int thread_id = 0) {

		if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);

		// TODO (lab 1) : if intersect with ray, use the returned information to compute the color ; otherwise black 
		// in lab 1, the color only includes direct lighting with shadows

		Vector P, N;
		double t;
		int object_id;
		if (intersect(ray, P, t, N, object_id)) {

			if (objects[object_id]->mirror) {
				Vector R = ray.u - 2.0 * dot(ray.u, N) * N;
				Ray R_ray(P + 1e-4 * N, R);
				// TO-DO: return getColor in the reflected direction, with recursion_depth+1 (recursively)
				return getColor(R_ray, recursion_depth + 1);
			} // else

			if (objects[object_id]->transparent) { // optional - fianlly implemented and it works!
				double n1 = 1.0;
				double n2 = 1.5;
				Vector N_real = N;

				if (dot(ray.u, N) > 0) {
					std::swap(n1, n2);
					N_real = -1.0 * N;
				}

				double dot_u_N = dot(ray.u, N_real);
				double cos_i = -dot_u_N;
				double rad = 1.0 - sqr(n1 / n2) * (1.0 - sqr(cos_i));
				Vector reflected_dir = ray.u - 2.0 * dot_u_N * N_real;
				if (rad < 0) {
					return getColor(Ray(P + 1e-4 * N_real, reflected_dir), recursion_depth + 1, thread_id);
				}
				double k0 = sqr((n1 - n2) / (n1 + n2));
				double R_fresnel = k0 + (1.0 - k0) * pow(1.0 - cos_i, 5.0);
				if (uniform(engine[thread_id]) < R_fresnel) {
					return getColor(Ray(P + 1e-4 * N_real, reflected_dir), recursion_depth + 1, thread_id);
				} else {
					Vector wTt = (n1 / n2) * (ray.u - dot_u_N * N_real);
					Vector wNt = -1.0 * N_real * sqrt(rad);
					Vector refracted_dir = wTt + wNt;
					return getColor(Ray(P - 1e-4 * N_real, refracted_dir), recursion_depth + 1, thread_id);
				}
			} // else

			// test if there is a shadow by sending a new ray
			// if there is no shadow, compute the formula with dot products etc.
			//CAREFUL.//Vector L = light_position - (1e-4 * N + P); //restore to original and ask about getting rid of the dsrk edges in mirror. 
			Vector L = light_position - P;  //Restored to original was trying to debug something.
			double d = L.norm();
			L.normalize();
			//Imp: Shadow ray, offset slightly to avoid self-intersection
			Ray shadow_ray(P + 1e-4 * N, L);
			Vector P_s, N_s;
			double t_s;
			int object_id_s;

			// If we hit an object and it's closer than the light (t_s < d), we are in shadow, so we need to implement that check here.
			double visibility = 1.0;
			if (intersect(shadow_ray, P_s, t_s, N_s, object_id_s) && t_s < d) {
				visibility = 0.0; // Pitch black shadow - new addition over the normal return method.
				//return Vector(0, 0, 0); // Pitch black shadow I assume?
			}
			// Compute diffuse shading (Lambertian) - Lab 2 --> include visibility float.
			double intensity = light_intensity / (4.0 * M_PI * d * d);
			double diffuse_factor = std::max(0.0, dot(N, L)); // Based on what we were taught in the lecture that we don't consider the negative stuff.
			Vector direct_lighting = objects[object_id]->albedo * (intensity * diffuse_factor * visibility / M_PI);

			// Lab 2 : Add indirect lighting component with a recursive Monte Carlo call
			// We need to randomly sample the hemisphere using a cosine-weighted distribution
			Vector indirect_direction = random_cos(N, thread_id);
			Ray randomRay(P + 1e-4 * N, indirect_direction);
			// Ray randomRay(P + 1e-4 * N, indirect_direction);
			// Only sum the indirect lighting if we haven't hit max recursion
			if (recursion_depth + 1 < max_light_bounce) {
				Vector indirect_lighting = objects[object_id]->albedo * getColor(randomRay, recursion_depth + 1, thread_id);
				return direct_lighting + indirect_lighting;
			}
			return direct_lighting;

			// Due to importance sampling, the Pi and Cos terms mathematically cancel out
			// Vector indirect_lighting = objects[object_id]->albedo * getColor(randomRay, recursion_depth + 1, thread_id);

			// // Note: Sir, I panicked a bit and couldn't recall this properly, according to me - Color = Intensity * Cos(theta) * (Albedo / Pi),  but as per notes eq: L = (I / 4*pi*d^2) * (rho / pi) * max(0, <N, w_i>), so I implemented what seemed more compliant to me.
			// Vector diffuse_color = objects[object_id]->albedo * (intensity * diffuse_factor / M_PI);

			//THIS WAS FOR LAB 1.
			// return diffuse_color;

			// TODO (lab 2) : add indirect lighting component with a recursive call
			//The methods defined above.
			//return direct_lighting + indirect_lighting; - Test out first.
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

	//I commented out the testing for the Lab 1 transparency part for now, since that is completely messign up.
	// Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), false);
	//Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), true, false); //check from the original.
	
	//Comment out the single thing.
	// Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8), true, false);

	// Create 3 spheres at different depths to demonstrate Depth of Field
	//Light intensity slightly darker because there are pure whites at many place - (darker or albedo smaller.)

	//Comment these out for Lab 3 or 4.

	// Sphere sphere_front(Vector(-12, 0, 25), 10., Vector(0.8, 0.3, 0.8), false); // Magenta - Blurry (Too close)
	// Sphere sphere_middle(Vector(0, 0, 0), 10., Vector(0.9, 0.9, 0.9), false, true);   // White - Sharp (At focal distance)
	// Sphere sphere_back(Vector(12, 0, -25), 10., Vector(0.2, 0.5, 0.8), true, false);  // Blue - Blurry (Too far)

	Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
	Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
	Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
	Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
	Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));

	//Be careful of the point of intersection and then it should normalize into the right thing when you are happy with it.
	TriangleMesh cat(Vector(0.8, 0.8, 0.8));
	cat.readOBJ("cat.obj");
    cat.scale_translate(0.6, Vector(0, -10, 0));
    cat.init_bvh();

	Scene scene;
	scene.camera_center = Vector(0, 0, 55); //Change in lecture.
	scene.light_position = Vector(-10,20,40);
	scene.light_intensity = 6E6; //3E7
	scene.fov = 60 * M_PI / 180.;  //I had previously set it to 55 * M_PI / 180.; too.
	scene.gamma = 2.2;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
	scene.max_light_bounce = 5;

	//Comment out the single sphere for now.
	//scene.addObject(&center_sphere);

	// scene.addObject(&sphere_front);
	// scene.addObject(&sphere_middle);
	// scene.addObject(&sphere_back);

	scene.addObject(&wall_left);
	scene.addObject(&wall_right);
	scene.addObject(&wall_front);
	scene.addObject(&wall_behind);
	scene.addObject(&ceiling);
	scene.addObject(&floor);
	
	//We msut add one for Mr. Cat too!
	scene.addObject(&cat);

	std::vector<unsigned char> image(W * H * 3, 0);


//
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; i++) {
		int thread_id = omp_get_thread_num(); // Thread safety for random generator engine
		for (int j = 0; j < W; j++) {

			Vector pixelColor(0., 0., 0.);
			int NB_PATHS = 64; // Defined as per the slides, effectively the number of rays we are sending per pixel, test more/less as I experimented later with 32, 64 and 1000 - multi-threaded hence it ran in 43 seconds.
			
			for (int k = 0; k < NB_PATHS; k++) {
				// Box-Muller Gaussian jitter per pixel for anti-aliasing
				double x_offset, y_offset;
				boxMuller(0.5, x_offset, y_offset, thread_id);
				
				// TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)
				double x = (j + x_offset) - W / 2.0 + 0.5;
				double y = H / 2.0 - (i + y_offset) - 0.5;
				double z = -W / (2.0 * tan(scene.fov / 2.0)); 
				
				Vector ray_direction(x, y, z);
				ray_direction.normalize();

				// Depth of field implementation - It was mentioned optional on the slides but I still wanted to implement it...
				double aperture_radius = 3.0; //Depth of Field - // The larger the radius, the stronger the DoF blur I assume...
				double focus_distance = 55.0; // Distance to the center_sphere (from Z=55 to Z=0....I f I remember correctly...)
				
				double t_focus = focus_distance / std::abs(ray_direction[2]);
				Vector focal_point = scene.camera_center + ray_direction * t_focus;
				
				double r1_lens = uniform(engine[thread_id]);
				double r2_lens = uniform(engine[thread_id]);
				double dx = aperture_radius * sqrt(r1_lens) * cos(2.0 * M_PI * r2_lens);
				double dy = aperture_radius * sqrt(r1_lens) * sin(2.0 * M_PI * r2_lens);
				
				Vector new_origin = scene.camera_center + Vector(dx, dy, 0);
				Vector new_direction = focal_point - new_origin;
				new_direction.normalize();

				//Commented out before Lab 3/4 - just for ease of notation.
				// Ray ray(new_origin, new_direction);

				//Old junk ray.
				// Ray ray(scene.camera_center, ray_direction);

				// pixelColor = pixelColor + scene.getColor(ray, 0, thread_id);
				pixelColor = pixelColor + scene.getColor(Ray(new_origin, new_direction), 0, thread_id);
			}

			// Average color from total paths
			pixelColor = pixelColor / NB_PATHS;

			image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(pixelColor[0] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(pixelColor[1] / 255., 1. / scene.gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(pixelColor[2] / 255., 1. / scene.gamma)));
		}
	}
	stbi_write_png("image_lab3_3.png", W, H, 3, &image[0], 0);

	return 0;
}

// Non-multi-threaded initial partial code - with the initial class intructions.
// #pragma omp parallel for schedule(dynamic, 1)
// 	for (int i = 0; i < H; i++) {
// 		int thread_id = omp_get_thread_num(); // I try to add thread safety for random generator engine, I hope it works.
// 		for (int j = 0; j < W; j++) {
// 			Vector color;

// 			// TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)			
// 			double x = j - W / 2.0 + 0.5; //Center horizontal pixel mapping.
// 			double y = H / 2.0 - i - 0.5; //Center the vertical mapping (basically invert the y-axis like explained int he class.)
// 			double z = -W / (2.0 * tan(scene.fov / 2.0)); //Map focal lenght via W and FOV.
// 			//New Ray and new origin remember.
			
// 			// Vector ray_direction(0., 0., -1); - Added the real x, y, z computed values now.
// 			Vector ray_direction(x, y, z);
// 			//Adding the ray direction normalization here too.
// 			ray_direction.normalize(); //we must ensure that the lenght is exactly 1 (reminder to self).


// 			Ray ray(scene.camera_center, ray_direction);

// 			// TODO (lab 2) : add Monte Carlo / averaging of random ray contributions here
// 			// TODO (lab 2) : add antialiasing by altering the ray_direction here
// 			// TODO (lab 2) : add depth of field effect by altering the ray origin (and direction) here

// 			color  = scene.getColor(ray, 0);

// 			image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(color[0] / 255., 1. / scene.gamma)));
// 			image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(color[1] / 255., 1. / scene.gamma)));
// 			image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(color[2] / 255., 1. / scene.gamma)));
// 		}
// 	}
// 	stbi_write_png("image.png", W, H, 3, &image[0], 0);

// 	return 0;
// }


//Notes for self during lecture and for future labs.
//Don't forget the object id and all too.
//Put white pixel if there is an intersection for example.
//don't forget to return the object_id.
//if image is totally white it might be possible your camera is inside the thing.
//you also wanna change the value and stuff to make sure things are considered.
//Don't forget to tweak the gamma function.
//Define your epsilon constant too, don't use lose values for those purposes its not the best programming practice.