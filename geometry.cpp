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

void print_color(color c){
    printf("%f, %f, %f\n", c.r, c.g, c.b);
}

void print_vector(point4d p){
    printf("%f, %f, %f, %f\n", p.x, p.y, p.z, p.w);
}

int quad_roots(GLdouble, GLdouble, GLdouble, GLdouble *);

const point4d operator*(const double lhs, const point4d &rhs){
    point4d p;
    p.x = lhs * rhs.x; 
    p.y = lhs * rhs.y;
    p.z = lhs * rhs.z;
    p.w = lhs * rhs.w;
    return p;
}
const point4d operator~(const point4d &rhs){
    point4d p;
    double mag = sqrt(rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z);
    p.x = rhs.x / mag;
    p.y = rhs.y / mag;
    p.z = rhs.z / mag;
    return p;
}

const color operator*(const double scale, const color &rhs){
    return {scale * rhs.r, scale * rhs.g, scale * rhs.b};
}

/*
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
            printf("Incorrect intersection: %f, %f, %f, %f, %f, %f, %f, "
                    "%f, %f, %f, %f, %f, %f, %f\n", a, b, c, d, e, f, g, h, j, k, zero, x, y, z);
            correct = false;


        }
    }
    return correct;
}
*/
#define randf(a) (((double)rand() - (double)RAND_MAX * 0.5) / (double)RAND_MAX * 2 * (double)a)

int herp(int argc, char **argv){
    vector4d u, v;
    u = {1, 6, 8, 0};
    v = {5, 4, 3, 0};
    vector4d reflect = u | (~v);
    print_vector(reflect);
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
        roots[0] = -1 * b / (2 * a);
        return 1;
    }else{ // 0 roots
        return 0;
    }
}

point4d* makepoint4d(GLdouble x, GLdouble y, GLdouble z) {
    point4d* p;
    /* allocate memory */
    p = (point4d*) malloc(sizeof(point4d));
    /* put stuff in it */
    p->x = x; p->y = y; p->z = z; 
    p->w = 1.0;
    return (p);
}

/* makes copy of point4d (or vector4d) */
point4d* copypoint4d(point4d *p0) {
    point4d* p;
    /* allocate memory */
    p = (point4d*) malloc(sizeof(point4d));

    p->x = p0->x;
    p->y = p0->y;
    p->z = p0->z;
    p->w = p0->w;         /* copies over vector4d or point4d status */
    return (p);
}

/* unallocates a point4d */
void freepoint4d(point4d *p) {
    if (p != NULL) {
        free(p);
    }
}

/* vector4d from point4d p to point4d q is returned in v */


void calculateDirection(point4d* p, point4d* q, point4d* v) {
    v->x = q->x - p->x;
    v->y = q->y - p->y;
    v->z = q->z - p->z;
    /* a direction is a point4d at infinity */
    v->w = 0.0;

    /* NOTE: v is not unit length currently, but probably should be */
}

/* given a vector4d, sets its contents to unit length */
void normalize(vector4d* v) {
    /* PUT YOUR CODE HERE */
}

/* point4d on ray r parameterized by t is returned in p */
void findpointOnRay(ray* r,double t,point4d* p) {
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
    s->c = makepoint4d(x,y,z);   /* center */
    s->r = r;   /* radius */
    s->m = NULL;   /* material */
    return(s);
}

/* returns TRUE if ray r hits sphere s, with parameter value in t */

int raySphereIntersect(ray* r,sphere* s,double* t) {
    point4d p;   /* start of transformed ray */
    double a,b,c;  /* coefficients of quadratic equation */
    double D;    /* discriminant */
    point4d* v;
    //p.x;
    /* transform ray so that sphere center is at origin */
    /* don't use matrix, just translate! */
    p.x = r->start->x - s->c->x;
    p.y = r->start->y - s->c->y;
    p.z = r->start->z - s->c->z;
    v = r->dir; /* point4d to direction vector4d */


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
        /* ignore roots which are less than zero (behind viewpoint4d) */
        if (*t < 0) {
            *t = (-b + D) / (2*a);
        }
        if (*t < 0) { return(FALSE); }
        else return(TRUE);
    }
}

/* normal vector4d of s at p is returned in n */
/* note: dividing by radius normalizes */
void findSphereNormal(sphere* s, point4d* p, vector4d* n) {
    n->x = (p->x - s->c->x) / s->r;  
    n->y = (p->y - s->c->y) / s->r;
    n->z = (p->z - s->c->z) / s->r;
    n->w = 0.0;
}


