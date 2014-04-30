/******************************************************************/
/*         Geometry functions                                     */
/*                                                                */
/* Group Members: <FILL IN>                                       */
/******************************************************************/

#ifdef _WIN32
#include <witdows.h>
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "raytrace.h"


int quad_roots(GLdouble, GLdouble, GLdouble, GLdouble *);

const point operator*(const double lhs, const point &rhs){
    point p;
    p.x = lhs * rhs.x; 
    p.y = lhs * rhs.y;
    p.z = lhs * rhs.z;
    return p;
}
const point operator~(const point &rhs){
    point p;
    double mag = sqrt(rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z);
    p.x = rhs.x / mag;
    p.y = rhs.y / mag;
    p.z = rhs.z / mag;
    return p;
}

class Object3D{
    private:
    public:
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
        //virtual int intersect(ray r, GLdouble intersects[4][4]) =0; // returns number of intersects
        virtual int intersect_t(ray r, GLdouble intersects[4]) =0; // returns parameter t for all intersections
};

class Polygon: public Object3D{
    public:
        // vertices on a face must be stored in counterclockwise order
        // otherwise the normals will be negated
        point *vertices;
        GLint num_vert;
        GLint **indices; // faces delimited by -1
        GLint *num_edges;
        GLint num_faces;
        Polygon(point *v, GLint nv, GLint **i, GLint *ne, GLint nf){
            vertices = v;
            num_vert = nv;
            indices = i;
            num_edges = ne;
            num_faces = nf;
        }
        /* copy constructor
           Polygon(const Polygon &other){
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
        int intersect_t(ray r, GLdouble intersects[4]){
            GLdouble enter = 0.0;
            GLdouble exit = 1000000.0;
            for(int i = 0; i < num_faces; i++){ // for every face
                vector v1 = vertices[indices[i][0]] - vertices[indices[i][1]];
                vector v2 = vertices[indices[i][2]] - vertices[indices[i][1]];
                vector normal = v1 ^ v2;
                GLdouble d = -1 * (normal * vertices[indices[i][0]]);
                // normal contains coefficients of x, y, z, and we have d, the constant
                // we have implicit form of plane, and 
                GLdouble face_dir = (*r.dir) * normal;
                if(face_dir > -.0001 && face_dir < 0.0001){
                    continue;
                }
                GLdouble t = -1 * (((*r.start) * normal) + d) / face_dir;
                if(face_dir < -0.0001){ // plane is facing towards us
                    if(t > enter)
                        enter = t;
                }else if(face_dir > 0.0001){
                    if(t < exit)
                        exit = t;
                }
            }
            if(enter < exit){
                intersects[0] = enter;
                intersects[1] = exit;
                return 2;
            }else if(enter == exit){
                intersects[0] = enter;
                return 1;
            }else
                return 0;
        }
};

class Quadric : public Object3D{ // class for general quadric solids
    public:
        GLdouble color[4];
        GLdouble a, b, c, d, e, f, g, h, j, k;
        Quadric(GLdouble a, GLdouble b, GLdouble c, GLdouble d, GLdouble e, GLdouble f, GLdouble g, GLdouble h, GLdouble j, GLdouble k, GLdouble color[4]){
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
            for(int i = 0; i < 4;i++){
                this -> color[i] = color[i];
            }
        }
        virtual int intersect_t(ray r, GLdouble t[4]){
            vector *v = r.dir;
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
        }
};

class Ellipsoid: public Quadric{
    public:
        Ellipsoid(GLdouble x, GLdouble y, GLdouble z, GLdouble xr, GLdouble yr, GLdouble zr, GLdouble color[4]): 
            Quadric(yr * zr, 0.0, 0.0,  xr * zr, 0.0, xr * yr,
                    -2 * x * yr * zr, -2 * y * xr * zr, -2 * z * xr * yr, x * x * yr * zr + y * y * xr * zr + z * z * xr * yr - xr * yr * zr, color){
                printf("yr * zr: %f\n", yr * zr);
            }
};

class Sphere : public Ellipsoid{
    public:
        Sphere(GLdouble x, GLdouble y, GLdouble z, GLdouble r2, GLdouble color[4]):
            Ellipsoid(x, y, z, r2, r2, r2, color){
                printf("x: %f\n", x);
            }
};

bool test_quadric_intersect(ray r, Quadric q){
    GLdouble intersects[4][4];
    int num_intersects = q.intersect(r, intersects);
    bool correct = true;
    GLdouble a = q.a;
    GLdouble b = q.b;
    GLdouble c = q.c;
    GLdouble d = q.d;
    GLdouble e = q.e;
    GLdouble f = q.f;
    GLdouble g = q.g;
    GLdouble h = q.h;
    GLdouble j = q.j;
    GLdouble k = q.k;
    GLdouble x, y, z, zero;
    for(int i = 0; i < num_intersects; i++){ 
        x = intersects[i][0];
        y = intersects[i][1];
        z = intersects[i][2];
        // plug intersect into quadric equation and check equality to 0
        if((zero = a * x * x + b * x * y + c * x * z + d * y * y + e * y * z +
                    f * z * z + g * x + h * y + j * z +
                    k) < -.0001 || zero > .0001){
            printf("Incorrect intersection: %f, %f, %f, %f, %f, %f, %f, \
                    %f, %f, %f, %f, %f, %f, %f\n", a, b, c, d, e, f, g, h, j, k, zero, x, y, z);
            correct = false;
        }
    }
    return correct;
}

#define randf(a) (((double)rand() - (double)RAND_MAX * 0.5) / (double)RAND_MAX * 2 * (double)a)
int main(int argc, char **argv){
    if(argc < 3){
        printf("Please enter seed and number of tests\n");
        return 1;
    }
    GLdouble color[4];
    srand(strtol(argv[1], NULL, 10));
    for(int i = 0; i < strtol(argv[2], NULL, 10); i++){
        Quadric q(randf(100), randf(100), randf(100), randf(100), randf(100), randf(100), randf(100), randf(100), randf(100), randf(100), color);
        point p(randf(100), randf(100), randf(100), 1.0);
        vector d(randf(100), randf(100), randf(100), 0.0);
        ray r(p, d); 
        test_quadric_intersect(r, q);
    }
    printf("Sphere\n");
    Sphere *sp = new Sphere(0.0, 0.0, 0.0, 64.0, color);
    Sphere s = *sp;
    printf("Sphere coefficients: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", s.a, s.b, s.c, s.d, s.e, s.f, s.g, s.h, s.j, s.k);
    printf("testing vector operators:\n");
    vector u(4, 5, 7, 0);
    vector v(-1, 6, 3, 0);
    printf("u: %f, %f, %f, %f ; v: %f, %f, %f, %f\n", u.x, u.y, u.z, u.w, v.x, v.y, v.z, v.w);
    vector result = u + v;
    printf("u + v: %f, %f, %f, %f\n", result.x, result.y, result.z, result.w);
    double dot = u * v;
    printf("u * v: %f\n", dot);
    result = u - v;
    printf("u - v: %f, %f, %f, %f\n", result.x, result.y, result.z, result.w);
    result = ~(u ^ v);
    printf("~(u ^ v): %f, %f, %f, %f\n", result.x, result.y, result.z, result.w);
    u -= v;
    printf("u -= v: %f, %f, %f, %f\n", u.x, u.y, u.z, u.w);
}
// Purely cubic solids will not be included. Since they can only
// have an odd number of real roots, there can exist no closed cubic
// surface.

// at^3 + bt^2 +ct + d = 0 can be tranformed to
// t^3 + pt + q = 0
// where p = (3 * a * c - b * b) / (3 * a * a)
// and q = (2 * pow(b, 3) - 9 * a * b * c + 27 * a * a * d) / (27 * pow(a, 3))

/*
   int cube_roots(GLdouble a, GLdouble b, GLdouble c, GLdouble d, GLdouble roots[3]){
   assert(a != 0.0);
   GLdouble p = (3 * a * c - b * b) / (3 * a * a);
   GLdouble q = (2 * pow(b, 3) - 9 * a * b * c + 27 * a * a * d)
   / (27 * pow(a, 3));
// t^3 + pt + q = 0, we now have this form
// substitute w - p / (3 * w) for t
// now we have w ^ 3 + q - p ^ 3 / (27 * w ^ 3);
// multiply by w ^ 3
// w ^ 6 + q * w ^ 3 - p ^ 3 / 27 = 0;


}
*/

// at^2 + bt^2 + c = 0
int quad_roots(GLdouble a, GLdouble b, GLdouble c, GLdouble roots[2]){
    // returns the number of roots 
    assert(a != 0.0); // make sure it is indeed a quadratic
    GLdouble discriminant = b * b - 4 * a * c;
    if(discriminant > 0.0001){ // 2 roots
        roots[0] = (-1 * b - sqrt(b * b - 4 * a * c)) / (2 * a);
        roots[1] = (-1 * b + sqrt(b * b - 4 * a * c)) / (2 * a);
        return 2;
    }else if(discriminant > - 0.0001 && discriminant < 0.0001){ // 1 root
        return 1;
    }else{ // 0 roots
        return 0;
    }
}

point* makePoint(GLdouble x, GLdouble y, GLdouble z) {
    point* p;
    /* allocate memory */
    p = (point*) malloc(sizeof(point));
    /* put stuff in it */
    p->x = x; p->y = y; p->z = z; 
    p->w = 1.0;
    return (p);
}

/* makes copy of point (or vector) */
point* copyPoint(point *p0) {
    point* p;
    /* allocate memory */
    p = (point*) malloc(sizeof(point));

    p->x = p0->x;
    p->y = p0->y;
    p->z = p0->z;
    p->w = p0->w;         /* copies over vector or point status */
    return (p);
}

/* unallocates a point */
void freePoint(point *p) {
    if (p != NULL) {
        free(p);
    }
}

/* vector from point p to point q is returned in v */
void calculateDirection(point* p, point* q, point* v) {
    v->x = q->x - p->x;
    v->y = q->y - p->y;
    v->z = q->z - p->z;
    /* a direction is a point at infinity */
    v->w = 0.0;

    /* NOTE: v is not unit length currently, but probably should be */
}

/* given a vector, sets its contents to unit length */
void normalize(vector* v) {
    /* PUT YOUR CODE HERE */
}

/* point on ray r parameterized by t is returned in p */
void findPointOnRay(ray* r,double t,point* p) {
    p->x = r->start->x + t * r->dir->x;
    p->y = r->start->y + t * r->dir->y;
    p->z = r->start->z + t * r->dir->z;
    p->w = 1.0;
}


/* SPHERES */

sphere* makeSphere(GLdouble x, GLdouble y, GLdouble z, GLdouble r) {
    sphere* s;
    /* allocate memory */
    s = (sphere*) malloc(sizeof(sphere));

    /* put stuff in it */
    s->c = makePoint(x,y,z);   /* center */
    s->r = r;   /* radius */
    s->m = NULL;   /* material */
    return(s);
}

/* returns TRUE if ray r hits sphere s, with parameter value in t */

int raySphereIntersect(ray* r,sphere* s,double* t) {
    point p;   /* start of transformed ray */
    double a,b,c;  /* coefficients of quadratic equation */
    double D;    /* discriminant */
    point* v;
    //p.x;
    /* transform ray so that sphere center is at origin */
    /* don't use matrix, just translate! */
    p.x = r->start->x - s->c->x;
    p.y = r->start->y - s->c->y;
    p.z = r->start->z - s->c->z;
    v = r->dir; /* point to direction vector */


    a = v->x * v->x  +  v->y * v->y  +  v->z * v->z;
    b = 2*( v->x * p.x  +  v->y * p.y  +  v->z * p.z);
    c = p.x * p.x + p.y * p.y + p.z * p.z - s->r * s->r;

    D = b * b - 4 * a * c;

    if (D < 0) {  /* no intersection */
        return (FALSE);
    }
    else {
        D = sqrt(D);
        /* First check the root with the lower value of t: */
        /* this one, since D is positive */
        *t = (-b - D) / (2*a);
        /* ignore roots which are less than zero (behind viewpoint) */
        if (*t < 0) {
            *t = (-b + D) / (2*a);
        }
        if (*t < 0) { return(FALSE); }
        else return(TRUE);
    }
}

/* normal vector of s at p is returned in n */
/* note: dividing by radius normalizes */
void findSphereNormal(sphere* s, point* p, vector* n) {
    n->x = (p->x - s->c->x) / s->r;  
    n->y = (p->y - s->c->y) / s->r;
    n->z = (p->z - s->c->z) / s->r;
    n->w = 0.0;
}


