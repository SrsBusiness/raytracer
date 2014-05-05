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
#include <stdarg.h>
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
color trace(ray*, Object3D *, int);
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
    ambi = {0.1, 0.1, 0.1};
    lsource = {1, 1, 1, 0};
    lsource = ~lsource;

    //s1 = makeSphere(0.0,0.0,-2.0,0.25);
    //s1->m = makeMaterial(1.0,0.1,0.15,0.3);
    //objects.push_back(new Sphere(0, 0, 2, 1, (material){1, 0, 0, .2}));
    objects.push_back(new Sphere(-.8, .2, -7, .25, (material){{1, 1, 1}, 0.2, 200.0, {0, .2, 0}, {0, 0, 1}}));
    objects.push_back(new Sphere(.6, .2, -5, .25, (material){{1, 1, 1}, 0.2, 200.0, {0, .2, 0}, {0, 1, 0}}));
    objects.push_back(new Sphere(0, .5, -8, .25, (material){{1, 1, 1}, 0.2, 200.0, {0, .2, 0}, {1, 0, 0}}));
    //objects.push_back(new Sphere(-.3, -.5, -4, .3, (material){{0.4, 0.3, 0.2}, {0, 0, 0}, 0.2, 5.0, {0, .2, 0}, {1, 0, 0}}));
    //objects.push_back(new Sphere(0, 0, -6, .5, (material){{0.4, 0.9, 0.1}, 0.1, 0.9, 2.0, 0.2, 0.4}));
    //objects.push_back(new Cube(0, 1, -6, .5, (material){{0.2, 0.7, 0.1}, 0.3, 0.6, 2.0, 0.2, 0.4}));

    //objects.push_back(new Cube(0, 0, 2, 1, (material){1, 0, 0, .7}));
}

void initCamera (int w, int h) {
    viewpoint = makepoint4d(0.0,0.0,0.0);
    pnear = 1.0;
    fovx = PI/6;
}

double max_d(int args, ...){
    va_list a;
    va_start(a, args);
    double max_d = va_arg(a, double);
    double next;
    for(int i = 1; i < args; i++){
        if((next = va_arg(a, double)) > max_d)
            max_d = next;
    }
    va_end(a);
    return max_d;
}

color cast(ray *);

void assert_color(color c){
    if(!(c.r >= 0.0 && c.g >= 0.0 && c.b >= 0.0)){
        printf("%f, %f, %f\n", c.r, c.g, c.b);
        return;
    }
    if(!(c.r <= 1.0 && c.g <= 1.0 && c.b <= 1.0)){
        printf("%f, %f, %f\n", c.r, c.g, c.b);
        return;
    }

}

void normalize(color *c){
    if(c -> r < 0.0)
        c -> r = 0.0;
    else if(c -> r > 1.0)
        c -> r = 1.0;
    if(c -> g < 0.0)
        c -> g = 0.0;
    else if(c -> g > 1.0)
        c -> g = 1.0;
    if(c -> b < 0.0)
        c -> b = 0.0;
    else if(c -> b > 1.0)
        c -> b = 1.0;
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
            //printf("primary ray\n");
            r.copy(worldPix, worldPix - *viewpoint);
            //printf("worldPix.w: %f, viewpoint.w: %f\n", worldPix.w, viewpoint -> w);
            //print_vector(*r.dir);
            eye = ~(*r.dir);

            /* trace the ray! */
            c = trace(&r, NULL, 0);
            //c = cast(&r);
            normalize(&c);
            //assert_color(c);
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
Object3D *firstHit(Object3D *exc, ray* r, point4d* p, vector4d* n, material *m) {
    GLdouble intersects[4];
    double t = 10000.0;
    point4d normals[4];
    point4d normal;
    Object3D *o = NULL; // object with which ray intersects first
    for(int i = 0; i < objects.size(); i++){
        if(exc != objects[i]){
            int num_intersects = objects[i] -> intersect_t(*r, intersects, normals);
            //printf("t: %f\n", t);
            if(num_intersects > 0 && intersects[0] < t && intersects[0] > 0.0){
                t = intersects[0];
                normal = normals[0];
                o = objects[i];
            }
        }
    }
    assert(t > 0.0);
    if(t < 10000.0){
        //printf("t: %f\n", t);
        findpointOnRay(r,t,p);
        //printf("point p: %f, %f, %f, %f\n", p -> x, p -> y, p -> z, p -> w);
        if (n != NULL) *n = normal;
        //printf("normal: %f, %f, %f, %f\n", normal.x, normal.y, normal.z, normal.w);
        if (m != NULL) *m = *(o -> m);
        return o;
    } else {
        p->w = 0.0;
        return NULL;
    }
}

color apply_material(GLdouble prop[3], color orig){
    return {prop[0] * orig.r, prop[1] * orig.g, prop[2] * orig.b};
}

color phong(point4d *p, vector4d *n, material *m, vector4d *in, vector4d *light, color *light_rgb){
    //printf("phong\n");
    vector4d l = ~(*light);
    vector4d normal = ~(*n);
    vector4d v = ~(-1 * *in);
    vector4d half = ~(l + v);
    vector4d reflected = ~(*in | normal);
    /* 
       printf("ambient\n");
       assert_color((m -> a * ambi));
       printf("diffuse\n");
       assert_color(((m -> d * max_d(2, 0.0, normal * l)) * source));
       printf("specular\n");
       assert_color(((m -> s * pow((normal * half) , m -> h)) * source));
       printf("m -> s: %f, pow((normal * half) , m -> h): %f\n", m -> s, pow((normal * half) , m -> h));
       printf("total\n");
       assert_color((m -> a * ambi) +
       ((m -> d * max_d(2, 0.0, normal * l)) * source) + 
       ((m -> s * pow((normal * half) , m -> h)) * source));
       */
    return m -> a * ambi +
        m -> d * (max_d(2, 0.0, normal * *light) * *light_rgb) +
        m -> s * (pow((max_d(2, 0.0, normal * half)) , m -> h) * *light_rgb);
}

color cast(ray *r){
    point4d intersect;
    vector4d normal;
    material m;
    firstHit(NULL, r, &intersect, &normal, &m);
    if(intersect.w == 0.0)
        return {0.0, 0.0, 0.0};
    //color result = phong(&intersect, &normal, &m, r -> dir, &lsource, &source);
    normal = ~normal;
    color result = {normal.x, normal.y, normal.z};
    //assert_color(result);
    return result;
}

color trace(ray* r, Object3D *prev, int depth) {
    ray flec, frac;
    color spec, refr, dull;
    color intensity;
    if(depth >= MAX_DEPTH){
        return back; 
    }
    //depth++;
    point4d intersect = { 0, 0, 0, 0 };
    vector4d n;
    material m;
    Object3D *obj = firstHit(prev, r, &intersect, &n, &m);
    n = ~n;
    if(intersect.w == 0){ // no intersect
        intensity = source * max_d(2, 0, (*(r->dir)) * lsource);
    } else {

        //printf("intersect\n");
        //print_vector(intersect);
        // compute reflection
        if (m.s.r + m.s.g + m.s.b > 0) {
            //printf("ray\n");
            //print_vector(~(*r -> dir));
            //printf("normal\n");
            //print_vector(n);
            vector4d norm = (*r->dir) | (~n);
            printf("depth: %d\n", depth);
            printf("Hit object\n");
            printf("normal\n");
            print_vector(n);
            print_color(m.d);
            //printf("reflect\n");
            //print_vector(norm);
            //print_vector((*r -> dir) | (~n));
            flec.copy(intersect, norm);
            //color incoming = trace(&flec, depth + 1);
            //spec = phong(&intersect, &n, &m, r -> dir, &norm, &incoming);
            color tmp = trace(&flec, obj, depth + 1);
            //vector4d half = ~(*r -> dir + n)
            //spec = apply_material(m.s, tmp) + source * pow((~(*r->dir) * ~eye), m.h);
            spec = m.s * tmp;
            //print_color(tmp);
            intensity = (spec + m.d) * .5;
            //spec = tmp * pow((~(*r->dir) * ~eye), m.h);
            /*
               printf("*r -> dir: ");
               print_vector(~*r -> dir);
               printf("eye: ");
               print_vector(~eye);
               printf("dot: %f\n", ~(*r->dir) * ~eye);
               */
            //print_color(source * pow((~(*r->dir) * ~eye), m.h));
            //
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
        //firstHit(&shadow, &shadow_sect, (vector4d*)nullptr, (material*)nullptr);
        /*
        if (shadow_sect.w == 0) {
            // no shadow (no intersection)
            //printf("hi\n");
            //dull = (source * apply_material(m.d * (n * lsource))) + (ambi * m.a);
        } else {
            //dull = ambi * m.a;
        }
        */
        dull = phong(&intersect, &n, &m, &eye, &lsource, &source); 
        intensity += dull;;
        /*
           if(intensity == spec){
           printf("spec\n");
           print_color(spec);
           printf("intensity before\n");
           print_color(intensity);
           }
           */
        //intensity = (intensity + m.d) * .5;
        //printf("material color\n");
        //print_color(m.d);
        //printf("intensity after\n");
        //print_color(intensity);
    }
    normalize(&intensity);
    return intensity;
}
