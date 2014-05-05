#ifndef _RAYTRACE_H_
#define _RAYTRACE_H_

/******************************************************************/
/*         Raytracer declarations                                 */
/******************************************************************/


/* constants */
#define TRUE 1
#define FALSE 0

#define PI 3.14159265358979323846264338327


/* data structures */

class Object3D;
class Polyhedron;
class Rect_prism;
class Cube;
class Quadric;
class Ellipsoid;
class Sphere;
class point4d;
class color;
const point4d operator*(const double, const point4d &);
const point4d operator~(const point4d &);
const color operator*(const double, const color &);


void print_vector(point4d);
void print_color(color);

class color {
    public:
        GLdouble r;
        GLdouble g;
        GLdouble b; 
        /* these should be between 0 and 1 */

        const color operator*(const GLdouble scale) const {
            return { r * scale, g * scale, b * scale };
        }

        const color operator*(const color& other){
            return {r * other.r, g * other.g, b * other.b};
        }

        const color operator*=(const GLdouble scale){
            *this = *this * scale;
            return *this;
        }

        const color operator+(const color &c) const {
            return { r + c.r, g + c.g, b + c.b };
        }

        const color operator+=(const color &c){
            *this = *this + c;
            return *this;
        }
        const bool operator==(const color &c){
            return r == c.r && g == c.g && b == c.b;
        }
        const bool operator!=(const color &c){
            return !(*this == c);
        }
};

typedef struct material {
    /* color */
    /* specular reflectivity */
    color s;
    /* refractivity */
    GLdouble r;
    /* shininess */
    GLdouble h;
    /* ambient reflectivity */
    color a;
    /* diffuse reflectivity */
    color d;
} material;

typedef struct sphere {
    point4d* c;  /* center */
    GLdouble r;  /* radius */
    material* m;
} sphere;


class point4d {
    public:
        GLdouble x;
        GLdouble y;
        GLdouble z;
        GLdouble w;
        point4d(GLdouble x, GLdouble y, GLdouble z, GLdouble w){
            this -> x = x;
            this -> y = y;
            this -> z = z;
            this -> w = w;
        }
        point4d(){
            x = y = z = w = 0;
        }
        // 2 * v

        const double operator*(const point4d &rhs) const { // dot product
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }
        const point4d operator*=(const double rhs){
            *this = rhs * *this; 
            return *this;
        }
        const point4d operator^(const point4d &rhs) const {
            point4d p;
            p.x = y * rhs.z - z * rhs.y;
            p.y = z * rhs.x - x * rhs.z;
            p.z = x * rhs.y - y * rhs.x;
            return p;
        }
        const point4d operator^=(const point4d &rhs){
            *this = *this ^ rhs;
            return *this;
        }
        const point4d operator+(const point4d &rhs) const {
            point4d p;
            p.x = x + rhs.x;
            p.y = y + rhs.y;
            p.z = z + rhs.z;
            p.w = w + rhs.w;
            return p;
        }
        const point4d operator +=(const point4d &rhs){
            *this = *this + rhs;
            return *this;
        }
        // r=d-2(d?n)n
        // reflection across normal
        const point4d operator|(const point4d &normal) const {
            //printf("operator reflect\n");
            //printf("this\n");
            //print_vector(*this);
            //printf("normal\n");
            //print_vector(normal);
            //printf("dot: %f\n", *this * normal);
            return *this - (2 *(*this * normal) * normal);
        }
        const point4d operator|=(const point4d &normal){
            *this = *this | normal;
            return *this;
        }
        const point4d operator-(const point4d &rhs) const {
            return (*this + (-1.0 * rhs)); 
        }
        const point4d operator-=(const point4d &rhs){
            *this = *this - rhs;
            return *this;
        }
        const point4d operator=(const point4d &rhs){
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            w = rhs.w;
            return *this;
        }
        const bool operator==(const point4d &rhs){
            return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w;
        }
        const bool operator!=(const point4d &rhs){
            return !(*this == rhs);
        }
};



/* a vector4d is just a point4d */
typedef point4d vector4d;

/* a ray is a start point4d and a direction */
typedef struct ray {
    point4d* start;
    vector4d* dir;
    ray(){
        start = new point4d(); 
        dir = new vector4d();
    }
    ray(point4d start, vector4d dir){
        this -> start = new point4d(start.x, start.y, start.z, start.w);
        this -> dir = new vector4d(dir.x, dir.y, dir.z, dir.w);
    }
    ray(const ray& other){
        this -> start = new point4d(other.start -> x, other.start -> y, other.start -> z, other.start -> w);
        this -> dir = new vector4d(other.dir -> x, other.dir -> y, other.dir -> z, other.dir -> w);
    }
    void copy(point4d start, vector4d dir){
        *this -> start = start;
        *this -> dir = dir;
    }
    ~ray(){
        delete start;
        delete dir;
    }
} ray;

int quad_roots(GLdouble, GLdouble, GLdouble, GLdouble[2]);

class Object3D{
    public:
        material *m;
        /*
           virtual int intersect(ray r, GLdouble intersects[4][4]){
        // substitute parametric equation of ray
        // into implicit quadric equation:
        // r: x = x0 + vxt
        //    y = y0 + vyt
        //    z = z0 + vzt
        GLdouble roots[4]; 
        GLdouble x = r.start -> x;
        GLdouble y = r.start -> y;
        GLdouble z = r.start -> z;
        GLdouble w = r.start -> w;
        int num_roots = intersect_t(r, roots);
        for(int i = 0; i < num_roots; i++){
        intersects[i][0] = x + r.dir -> x * roots[i];
        intersects[i][1] = y + r.dir -> y * roots[i];
        intersects[i][2] = z + r.dir -> z * roots[i];
        intersects[i][3] = w + r.dir -> w * roots[i]; // r.dir -> w should be 0, and w will remain 1
        }
        return num_roots;
        }
        */
        //virtual int intersect(ray r, GLdouble intersects[4][4]) =0; // returns number of intersects
        virtual int intersect_t(ray, GLdouble[4], point4d[4]) =0; // returns parameter t for all intersections
};

class Polyhedron: public Object3D{
    public:
        // vertices on a face must be stored in counterclockwise order
        // otherwise the normals will be negated
        point4d *vertices;
        GLint num_vert;
        GLint **indices; // faces delimited by -1
        GLint *num_edges;
        GLint num_faces;
        Polyhedron(point4d *v, GLint nv, GLint **ind, GLint *ne, GLint nf, material m){
            num_vert = nv;
            vertices = new point4d[num_vert];
            for(int i = 0; i < num_vert; i++)
                vertices[i] = v[i];
            num_faces = nf;
            num_edges = new GLint[num_faces];
            indices = new GLint *[num_faces];
            for(int i = 0; i < num_faces; i++){
                num_edges[i] = ne[i];
                indices[i] = new GLint[num_edges[i]];
                for(int j = 0; j < num_edges[i]; j++)
                    indices[i][j] = ind[i][j];
            }
            this -> m = new material();
            *this -> m = m;
        }
        Polyhedron(GLint nv, GLint nf, material m){
            num_vert = nv;
            num_faces = nf;
            vertices = new point4d[num_vert];
            num_edges = new GLint[num_faces];
            indices = new GLint *[num_faces];
            this -> m = new material();
            *this -> m = m;
        }
        ~Polyhedron(){
            delete vertices;
            delete num_edges;
            for(int i = 0; i < num_faces; i++){
                delete indices[i];
            }
            delete indices;
        }
        /* copy constructor
           Polyhedron(const Polyhedron &other){
           num_vert = other.num_vert;
           vertices = new GLdouble[num_vert];
           for(int i = 0; i < num_vert; i++)
           vertices[i] = other.vertices[i];
           num_faces = other.num_faces;
           index_len  = other.index_len;
           indices = new GLint[index_len];
           for(int i = 0; i < index_len; i++)
           indices[i] = other.indices[i]; 
           }
           */
        int intersect_t(ray r, GLdouble intersects[4], point4d normals[4]){
            GLdouble enter = 0.0;
            GLdouble exit = 1000000.0;
            point4d enter_n, exit_n;
            for(int i = 0; i < num_faces; i++){ // for every face
                vector4d v1 = vertices[indices[i][0]] - vertices[indices[i][1]];
                vector4d v2 = vertices[indices[i][2]] - vertices[indices[i][1]];
                vector4d normal = v2 ^ v1;
                GLdouble d = -1 * (normal * vertices[indices[i][0]]);
                // normal contains coefficients of x, y, z, and we have d, the constant
                // we have implicit form of plane, and 
                GLdouble face_dir = (*r.dir) * normal;
                if(face_dir > -.0001 && face_dir < 0.0001){
                    continue;
                }
                GLdouble t = -1 * (((*r.start) * normal) + d) / face_dir;
                if(face_dir < -0.0001){ // plane is facing towards us
                    if(t > enter){
                        enter = t;
                        enter_n = normal;
                    }
                }else if(face_dir > 0.0001){
                    if(t < exit)
                        exit = t;
                    exit_n = normal;
                }
            }
            if(enter < exit){
                intersects[0] = enter;
                intersects[1] = exit;
                normals[0] = enter_n;
                normals[1] = exit_n;
                return 2;
            }else if(enter == exit){
                intersects[0] = enter;
                normals[0] = enter_n;
                normals[1] = exit_n;
                return 1;
            }else
                return 0;
        }
};

class Rect_Prism : public Polyhedron{
    public:
        Rect_Prism(GLdouble x, GLdouble y, GLdouble z, GLdouble w, GLdouble h, GLdouble l, material m) : Polyhedron(8, 6, m){
            vertices[0] = {x, y, z, 0};
            vertices[1] = {x + w, y, z, 0};
            vertices[2] = {x + w, y, z + l, 0};
            vertices[3] = {x, y, z + l, 0};
            vertices[4] = {x, y + h, z, 0};
            vertices[5] = {x + w, y + h, z, 0};
            vertices[6] = {x + w, y + h, z + l, 0};
            vertices[7] = {x, y + h, z + l, 0};
            for(int i = 0; i < num_faces; i++)
                num_edges[i] = 4;
            for(int i = 0; i < num_faces; i++){
                indices[i] = new int[num_edges[i]];
            }
            indices[0][0] = 0; indices[0][1] = 1; indices[0][2] = 2; indices[0][3] = 3;
            indices[1][0] = 1; indices[1][1] = 5; indices[1][2] = 6; indices[1][3] = 2;
            indices[2][0] = 0; indices[2][1] = 4; indices[2][2] = 5; indices[2][3] = 1;
            indices[3][0] = 3; indices[3][1] = 7; indices[3][2] = 4; indices[3][3] = 0;
            indices[4][0] = 2; indices[4][1] = 6; indices[4][2] = 7; indices[4][3] = 3;
            indices[5][0] = 6; indices[5][1] = 5; indices[5][2] = 4; indices[5][3] = 7;
        }
};

class Cube: public Rect_Prism{
    public:
        Cube(GLdouble x, GLdouble y, GLdouble z, GLdouble side, material m) : Rect_Prism(x, y, z, side, side, side, m){
        }
};

class Quadric : public Object3D{ // class for general quadric solids
    public:
        GLdouble a, b, c, d, e, f, g, h, j, k;
        Quadric(GLdouble a, GLdouble b, GLdouble c, GLdouble d, GLdouble e, GLdouble f, GLdouble g, GLdouble h, GLdouble j, GLdouble k, material m){
            this -> a = a;
            this -> b = b;
            this -> c = c;
            this -> d = d;
            this -> e = e;
            this -> f = f;
            this -> g = g;
            this -> h = h;
            this -> j = j;
            this -> k = k;
            this -> m = new material();
            *this -> m = m;
        }
        virtual int intersect_t(ray r, GLdouble t[4], point4d normals[4]){
            vector4d *v = r.dir;
            GLdouble x = r.start -> x;
            GLdouble y = r.start -> y;
            GLdouble z = r.start -> z;
            GLdouble w = r.start -> w;
            GLdouble a_q = a * v -> x * v -> x + b * v -> x * v -> y +
                c * v -> x * v -> z + d * v -> y * v -> y +
                e * v -> y * v -> z + f * v -> z * v -> z;
            GLdouble b_q =  2 * a * x * v -> x + b * (x * v -> y + y * v -> x) +
                c * (x * v -> z + z * v -> x) + 2 * d * y * v -> y +
                e * (y * v -> z + z * v -> y) + 2 * f * z * v -> z +
                g * v -> x + h * v -> y + j * v -> z;
            GLdouble c_q = a * x * x + b * x * y + c * x * z + d * y * y +
                e * y * z + f * z * z + g * x + h * y + j * z + k;
            //GLdouble roots[2]; // solutions for t
            int num_roots = quad_roots(a_q, b_q, c_q, t);
            // normals
            // n.x = 2ax + by + cz + g
            // n.y = 2dy + bx + ez + h
            // n.z = 2fz + cx + ey + j
            for(int i = 0; i < num_roots; i++){
                GLdouble tmp_x = x + v -> x * t[i];
                GLdouble tmp_y = y + v -> y * t[i];
                GLdouble tmp_z = z + v -> z * t[i];
                normals[i] = {2 * a * tmp_x + b * tmp_y + c * tmp_z + g, 
                    2 * d * tmp_y + b * tmp_x + e * tmp_z + h, 
                    2 * f * tmp_z + c * tmp_x + e * tmp_y + j,
                    0};
            }
            return num_roots;
        }
};

class Ellipsoid: public Quadric{
    public:
        Ellipsoid(GLdouble x, GLdouble y, GLdouble z, GLdouble xr, GLdouble yr, GLdouble zr, material m): 
            Quadric(yr * zr, 0.0, 0.0,  xr * zr, 0.0, xr * yr,
                    -2 * x * yr * zr, -2 * y * xr * zr, -2 * z * xr * yr, x * x * yr * zr + y * y * xr * zr + z * z * xr * yr - xr * yr * zr, m){
            }
};

class Sphere : public Ellipsoid{
    public:
        Sphere(GLdouble x, GLdouble y, GLdouble z, GLdouble r2, material m):
            Ellipsoid(x, y, z, r2, r2, r2, m){
            }
};



/* functions in raytrace.cpp */
void traceRay(ray*, color*, int);

/* functions in geometry.cpp */
sphere* makeSphere(GLdouble, GLdouble, GLdouble, GLdouble);
point4d* makepoint4d(GLdouble, GLdouble, GLdouble);
point4d* copypoint4d(point4d *);
void freepoint4d(point4d *);
void calculateDirection(point4d*,point4d*,point4d*);
void findpointOnRay(ray*,double,point4d*);
int raySphereIntersect(ray*,sphere*,double*);
void findSphereNormal(sphere*,point4d*,vector4d*);

/* functions in light.cpp */
material* makeMaterial(GLdouble, GLdouble, GLdouble, GLdouble);
void shade(point4d*,vector4d*,material*,vector4d*,color*,int);

/* global variables */
extern int width;
extern int height;

#endif	/* _RAYTRACE_H_ */
