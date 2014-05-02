/******************************************************************/
/*         Main raytracer file                                    */
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
#include "lowlevel.h"
#include "raytrace.h"
#include <vector>

#define MAX_DEPTH 5
/* local functions */
void initScene(void);
void initCamera (int, int);
void display(void);
void init(int, int);
void traceRay(ray*,color*);
void drawScene(void);
void firstHit(ray*,point4d*,vector4d*,material**);
color trace(ray*,int);
color normalized(color);

using namespace std;
/* local data */

color source, back, ambi;
vector4d lsource;
vector4d eye;

/* the scene: so far, just one sphere */
//sphere* s1;
vector<Object3D *> objects;
/* the viewing parameters: */
point4d* viewpoint;
GLdouble pnear;  /* distance from viewpoint4d to image plane */
GLdouble fovx;  /* x-angle of view frustum */
int width = 500;     /* width of window in pixels */
int height = 350;    /* height of window in pixels */

int main (int argc, char** argv) {
    int win;

    glutInit(&argc,argv);
    glutInitWindowSize(width,height);
    glutInitWindowPosition(100,100);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    win = glutCreateWindow("raytrace");
    glutSetWindow(win);
    init(width,height);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}

void init(int w, int h) {

    /* OpenGL setup */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
    glClearColor(0.0, 0.0, 0.0, 0.0);  

    /* low-level graphics setup */
    initCanvas(w,h);

    /* raytracer setup */
    initCamera(w,h);
    initScene();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    drawScene();  /* draws the picture in the canvas */
    flushCanvas();  /* draw the canvas to the OpenGL window */
    glFlush();
}

void initScene () {
    source = {1, 1, 1};
    back = {0, 0, 0};
    ambi = {0.2, 0.1, 0.1};
    lsource = {0.5, 0.25, 0.75, 0};
    lsource = ~lsource;

    //s1 = makeSphere(0.0,0.0,-2.0,0.25);
    //s1->m = makeMaterial(1.0,0.1,0.15,0.3);
    //objects.push_back(new Sphere(0, 0, 2, 1, (material){1, 0, 0, .2}));
    objects.push_back(new Sphere(0, 0, -5, .25, (material){{0.4, 0.3, 0.2}, 0.5, 0.2, 0.5, 0.2, 0.4}));
    objects.push_back(new Sphere(.1, 0, -4.7, .25, (material){{0.4, 0.9, 0.1}, 0.1, 0.9, 0.8, 0.2, 0.4}));
    objects.push_back(new Cube(0, 1, -6, .5, (material){{0.2, 0.7, 0.1}, 0.3, 0.6, 0.2, 0.2, 0.4}));

    //objects.push_back(new Cube(0, 0, 2, 1, (material){1, 0, 0, .7}));
}

void initCamera (int w, int h) {
    viewpoint = makepoint4d(0.0,0.0,0.0);
    pnear = 1.0;
    fovx = PI/6;
}

void drawScene () {
    int i, j;
    GLdouble imageWidth;
    /* declare data structures on stack to avoid dynamic allocation */
    point4d worldPix;  /* current pixel in world coordinates */
    point4d direction; 
    ray r;
    color c;

    /* initialize */
    worldPix.w = 1.0;
    worldPix.z = -pnear;

    imageWidth = 2 * pnear * tan(fovx / 2);

    /* trace a ray for every pixel */
    for (i = 0; i < width; i++) {
       /* Refresh the display */
       /* Comment this line out after debugging */
       flushCanvas();

       for (j = 0; j < height; j++) {

           /* find position of pixel in world coordinates */
           /* y position = (pixel height/middle) scaled to world coords */ 
           worldPix.y = (j - (height / 2)) * imageWidth / width;
           /* x position = (pixel width/middle) scaled to world coords */ 
           worldPix.x = (i - (width / 2)) * imageWidth / width;

           /* find direction */
           /* note: direction vector4d is NOT NORMALIZED */
           r.copy(worldPix, worldPix - *viewpoint);
           eye = ~(*r.dir);

           /* trace the ray! */
           c = trace(&r, 0);
           /* write the pixel! */
           drawPixel(i,j,c.r,c.g,c.b);
       }
    }
}

/* returns the color seen by ray r in parameter c */
/* d is the recursive depth */

//void traceRay(ray* r, color* c, int d) {
//    point4d p;  /* first intersection point4d */
//    vector4d n;
//    material* m;
//
//    p.w = 0.0;  /* inialize to "no intersection" */
//    firstHit(r,&p,&n,&m);
//    if (p.w != 0.0) {
//        shade(&p, &n, m, r->dir, c, d);  /* do the lighting calculations */
//    } else {             /* nothing was hit */
//        //printf("hit nothing\n");
//        c->r = 0.0;
//        c->g = 0.0;
//        c->b = 0.0;
//    }
//}



/* firstHit */
/* If something is hit, returns the finite intersection point4d p, 
   the normal vector4d n to the surface at that point4d, and the surface
   material m. If no hit, returns an infinite point4d (p->w = 0.0) */
void firstHit(ray* r, point4d* p, vector4d* n, material *m) {
    GLdouble intersects[4];
    double t = 10000.0;
    point4d normals[4];
    point4d normal;
    Object3D *o = NULL; // object with which ray intersects first
    for(int i = 0; i < objects.size(); i++){
        int num_intersects = objects[i] -> intersect_t(*r, intersects, normals);
        //printf("t: %f\n", t);
        if(num_intersects > 0 && intersects[0] < t){
            t = intersects[0];
            normal = normals[0];
            o = objects[i];
        }else{
            //printf("no intersect\n");
        }
    }
    if(t < 10000.0){
        //printf("t: %f\n", t);
        findpointOnRay(r,t,p);
        //printf("point p: %f, %f, %f, %f\n", p -> x, p -> y, p -> z, p -> w);
        if (n != NULL) *n = normal;
        //printf("normal: %f, %f, %f, %f\n", normal.x, normal.y, normal.z, normal.w);
        if (m != NULL) *m = *(o -> m);
    } else {
        p->w = 0.0;
    }
}

color trace(ray* r, int depth) {
    ray flec, frac;
    color spec, refr, dull;
    color intensity;
    if(depth >= MAX_DEPTH){
        return back; 
    }
    depth++;
    point4d intersect = { 0, 0, 0, 0 };
    vector4d n;
    material m;
    firstHit(r, &intersect, &n, &m);
    if(intersect.w == 0){ // no intersect
        intensity = source * (*(r->dir) * lsource);
    } else {
        // compute reflection
        if (m.s > 0) {
            vector4d norm = (*r->dir)|n;
            flec.copy(intersect, norm); 
            spec = trace(&flec, depth) * m.s + source * pow((~(*r->dir) * eye), m.h);
        } else {
            spec = { 0, 0, 0 };
        }
        
        // compute refraction
        if (m.r > 0) {
            refr = { 0, 0, 0 };
        } else {
            refr = { 0, 0, 0 };
        }

        // compute shadow
        ray shadow;
        shadow.copy(intersect, lsource);
        point4d shadow_sect = { 0, 0, 0, 0 };
        firstHit(&shadow, &shadow_sect, (vector4d*)nullptr, (material*)nullptr);
        if (shadow_sect.w == 0) {
            // no shadow (no intersection)
            dull = (source * (m.d * (n * lsource))) + (ambi * m.a);
        } else {
            dull = ambi * m.a;
        }

        intensity = spec + refr + dull;
    }
    return intensity;
}
