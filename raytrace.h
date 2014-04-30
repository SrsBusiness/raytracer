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

class object;
class point;
const point operator*(const double, const point &);
const point operator~(const point &);

class point {
    public:
        GLdouble x;
        GLdouble y;
        GLdouble z;
        GLdouble w;
        point(GLdouble x, GLdouble y, GLdouble z, GLdouble w){
            this -> x = x;
            this -> y = y;
            this -> z = z;
            this -> w = w;
        }
        point(){
            x = y = z = w = 0;
        }
        // 2 * v

        const double operator*(const point &rhs){ // dot product
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }
        const point operator*=(const double rhs){
            *this = rhs * *this; 
            return *this;
        }
        const point operator^(const point &rhs){
            point p;
            p.x = y * rhs.z - z * rhs.y;
            p.y = z * rhs.x - x * rhs.z;
            p.z = x * rhs.y - y * rhs.x;
            return p;
        }
        const point operator^=(const point &rhs){
            *this = *this ^ rhs;
            return *this;
        }
        const point operator+(const point &rhs){
            point p;
            p.x = x + rhs.x;
            p.y = y + rhs.y;
            p.z = z + rhs.z;
            p.w = w + rhs.w;
            return p;
        }
        const point operator +=(const point &rhs){
            *this = *this + rhs;
            return *this;
        }
        const point operator-(const point &rhs){
            return (*this + (-1.0 * rhs)); 
        }
        const point operator-=(const point &rhs){
            *this = *this - rhs;
            return *this;
        }
        
};



/* a vector is just a point */
typedef point vector;

/* a ray is a start point and a direction */
typedef struct ray {
    point* start;
    vector* dir;
    ray(){
    }
    ray(point start, vector dir){
        this -> start = new point(start.x, start.y, start.z, start.w);
        this -> dir = new vector(dir.x, dir.y, dir.z, dir.w);
    }
    ray(const ray& other){
        this -> start = new point(other.start -> x, other.start -> y, other.start -> z, other.start -> w);
        this -> dir = new vector(other.dir -> x, other.dir -> y, other.dir -> z, other.dir -> w);
    }
    ~ray(){
        delete start;
        delete dir;
    }
} ray;

typedef struct material {
    /* color */
    GLdouble r;
    GLdouble g;
    GLdouble b; 
    /* ambient reflectivity */
    GLdouble amb;
} material;

typedef struct color {
    GLdouble r;
    GLdouble g;
    GLdouble b; 
    /* these should be between 0 and 1 */
} color;

typedef struct sphere {
    point* c;  /* center */
    GLdouble r;  /* radius */
    material* m;
} sphere;

/* functions in raytrace.cpp */
void traceRay(ray*, color*, int);

/* functions in geometry.cpp */
sphere* makeSphere(GLdouble, GLdouble, GLdouble, GLdouble);
point* makePoint(GLdouble, GLdouble, GLdouble);
point* copyPoint(point *);
void freePoint(point *);
void calculateDirection(point*,point*,point*);
void findPointOnRay(ray*,double,point*);
int raySphereIntersect(ray*,sphere*,double*);
void findSphereNormal(sphere*,point*,vector*);

/* functions in light.cpp */
material* makeMaterial(GLdouble, GLdouble, GLdouble, GLdouble);
void shade(point*,vector*,material*,vector*,color*,int);

/* global variables */
extern int width;
extern int height;

#endif	/* _RAYTRACE_H_ */
