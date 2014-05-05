/******************************************************************/
/*         Lighting functions                                     */
/*                                                                */
/* Group Members: <FILL IN>                                       */
/******************************************************************/

#ifdef _WIN32
#include <windows.h>
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "raytrace.h"
/*
material* makeMaterial(GLdouble r, GLdouble g, GLdouble b, GLdouble amb) {
    material* m;

    m = (material*) malloc(sizeof(material));
    m->c = { r, g, b };
    m->a = amb;
    return(m);
}
*/

/* LIGHTING CALCULATIONS */

/* shade */
/* color of point4d p with normal vector4d n and material m returned in c */
/* in is the direction of the incoming ray and d is the recusive depth */
/*
void shade(point4d* p, vector4d* n, material* m, vector4d* in, color* c, int d) {

    c->r = m->a * m->c.r;
    c->g = m->a * m->c.g;
    c->b = m->a * m->c.b;

    if (c->r > 1.0) c->r = 1.0;
    if (c->g > 1.0) c->g = 1.0;
    if (c->b > 1.0) c->b = 1.0;

}
*/
