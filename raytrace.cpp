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
#include <pthread.h>

#define MAX_DEPTH 16
#define E 2.71828182846
#define IN 0
#define OUT 1

#define randfs(a) (((double)rand() - (double)RAND_MAX * 0.5) / (double)RAND_MAX * 2 * (double)a) // signed
#define randf(a) ((double)rand() / (double)RAND_MAX * (double)a)
//#define MAX_THREADS 8
/* local functions */
void initScene(void);
void initCamera (int, int);
void display(void);
void init(int, int);
void traceRay(ray*,color*);
void drawScene(void);
Object3D *firstHit(Object3D *exc, ray* r, point4d* p, vector4d* n, material *m, int);
double trace(ray*, Object3D *, int, int, int);
color normalized(color);

using namespace std;
/* local data */

color source, back, ambi;
vector4d lsource;

bool traced = false;

int max_threads = 1;
int samples = 1;
/* the scene: so far, just one sphere */
//sphere* s1;
vector<Object3D *> objects;
/* the viewing parameters: */
point4d* viewpoint;
GLdouble pnear;  /* distance from viewpoint4d to image plane */
GLdouble fovx;  /* x-angle of view frustum */
int width = 1200;     /* width of window in pixels */
int height = 800;    /* height of window in pixels */
GLdouble w_length[3] = {RED * RED, GREEN * GREEN, BLUE * BLUE};
double gaussian(double, double);

int main (int argc, char** argv) {
    /*
       double d = sqrt(2);
       double x = gaussian(.30, .40);
       printf("gaussian(%f, %f) = %f\n", .30, .40, x);
       x = gaussian(.40, .30);
       printf("gaussian(%f, %f) = %f\n", .40, .30, x);
       x = gaussian(.50, .0);
       printf("gaussian(%f, %f) = %f\n", .50, .00, x);
       for(int i = 0; i < 1000; i++){
       printf("%f\n", randfs(.5));
       }
       */
    if(argc < 2){
        printf("Please specify number of threads\n");
        return 1;
    }
    max_threads = strtol(argv[1], NULL, 10);
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
    back = {0.0, 0.0, 0.0};
    ambi = {0.0, 0.0, 0.0};
    lsource = {0, 0, 1, 0};
    lsource = ~lsource;
    double r, g, b;
    objects.push_back(new Sphere(-1, -1, -11, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.0, 0.0, 0.8}, {0.0, 0.0, 0.2}}));
    objects.push_back(new Sphere(-1, 0, -12, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.0, 0.8, 0.0}, {0.0, 0.2, 0.0}}));
    objects.push_back(new Sphere(-1, 1, -11, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.0, 0.8, 0.8}, {0.0, 0.2, 0.2}}));
    
    objects.push_back(new Sphere(0, -1, -12, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.8, 0.0, 0.0}, {0.2, 0.0, 0.0}}));
    objects.push_back(new Sphere(0, 0, -11, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.8, 0.0, 0.8}, {0.2, 0.0, 0.3}}));
    objects.push_back(new Sphere(0, 1, -12, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.8, 0.8, 0.0}, {0.2, 0.2, 0.0}}));
    
    objects.push_back(new Sphere(1, -1, -11, .25, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.8, 0.8, 0.8}, {0.2, 0.2, 0.2}}));
    objects.push_back(new Sphere(0, 0, -10, .1, (material){{1, 1, 1}, 1.4580, .00354, 1000.0, {0.8, 0.8, 0.8}, {0.2, 0.2, 0.2}}));
    /*
    for(int i = -2; i <= 2; i++){
        for(int j = -1; j <= 1; j++){
            r = randf(1.0);
            g = randf(1.0);
            b = randf(1.0);
            double z;
            if((j + i) & 1)
                z = -11.0;
            else
                z = -12.0;
            objects.push_back(new Sphere(i, j, z, .25, (material){{1, 1, 1}, 0.2, 1000.0, {.5, .5, .5}, {.5, .5, .5}}));
        }
    }
    for(double i = -1.0; i <= 1.0; i += .5){
        for(double j = -1.0; j <= 1.0; j += .5){
            r = randf(1.0);
            g = randf(1.0);
            b = randf(1.0);
            objects.push_back(new Sphere(i + .25, j + .25, -8, 0.05, (material){{1, 1, 1}, 0.2, 1000.0, {.5, .5, .5}, {.5, .5, .5}}));
        }
    }
    */
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
    if(!(c.rgb[0] >= 0.0 && c.rgb[1] >= 0.0 && c.rgb[2] >= 0.0)){
        printf("%f, %f, %f\n", c.rgb[0], c.rgb[1], c.rgb[2]);
        return;
    }
    if(!(c.rgb[0] <= 1.0 && c.rgb[1] <= 1.0 && c.rgb[2] <= 1.0)){
        printf("%f, %f, %f\n", c.rgb[0], c.rgb[1], c.rgb[2]);
        return;
    }

}

void normalize(color *c){
    if(c -> rgb[0] < 0.0)
        c -> rgb[0] = 0.0;
    else if(c -> rgb[0] > 1.0)
        c -> rgb[0] = 1.0;
    if(c -> rgb[1] < 0.0)
        c -> rgb[1] = 0.0;
    else if(c -> rgb[1] > 1.0)
        c -> rgb[1] = 1.0;
    if(c -> rgb[2] < 0.0)
        c -> rgb[2] = 0.0;
    else if(c -> rgb[2] > 1.0)
        c -> rgb[2] = 1.0;
}

GLdouble gaussian(GLdouble x, GLdouble y){
    return pow(E, (x * x + y * y) / -2) / (2 * PI);
}

void *draw(void *ptr){
    long thread_id = (long)ptr;
    GLdouble imageWidth;
    point4d worldPix;  /* current pixel in world coordinates */
    point4d direction; 
    ray r;
    color c;

    worldPix.w = 1.0;
    worldPix.z = -pnear;

    double total;
    imageWidth = 2 * pnear * tan(fovx / 2);
    for(int i = thread_id * width / max_threads; i < width && i < (thread_id + 1) * width / max_threads; i++){
        for(int j = 0; j < height; j++){
            //total = 0.0;
            worldPix.y = (j - (height / 2)) * imageWidth / width;
            worldPix.x = (i - (width / 2)) * imageWidth / width;
            r.copy(worldPix, worldPix - *viewpoint);
            double gauss = gaussian(0.0, 0.0);
            c = (color){{trace(&r, NULL, 0, 0, OUT), trace(&r, NULL, 0, 1, OUT), trace(&r, NULL, 0, 2, OUT)} } * gauss;
            total = gauss;
            //printf("samples: \n");
            for(int k = 0; k < samples; k++){ 
                double dx = randfs(.5);
                double dy = randfs(.5);
                gauss = gaussian(dx, dy);
                point4d sample;
                sample.y = ((double)j + dy - (height / 2)) * imageWidth / width;
                sample.x = ((double)i + dx - (width / 2)) * imageWidth / width;
                sample.z = -pnear;
                sample.w = 1.0;
                //printf("worldPix.x: %f, worldPix.y: %f, sample.x: %f, sample.y: %f\n", worldPix.x, worldPix.y, sample.x, sample.y);
                r.copy(sample, sample - *viewpoint);
                c += (color){{trace(&r, NULL, 0, 0, OUT), trace(&r, NULL, 0, 1, OUT), trace(&r, NULL, 0, 2, OUT)}} * gauss;
                total += gauss;
                //normalize(&c);
            }
            c /= total;
            drawPixel(i,j,c.rgb[0],c.rgb[1],c.rgb[2]);
        }
    }
}

void drawScene () {
    if(!traced){
        pthread_t threads[max_threads];
        for(long i = 0; i < max_threads; i++){
            pthread_create(threads + i, NULL, draw, (void *)i);
        }
        /* trace a ray for every pixel */
        for(int i = 0; i < max_threads; i++){
            pthread_join(threads[i], NULL);
        }
        traced = true;
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
Object3D *firstHit(Object3D *exc, ray* r, point4d* p, vector4d* n, material *m, int in_out) {
    GLdouble intersects[4];
    double t = 10000.0;
    point4d normals[4];
    point4d normal;
    Object3D *o = NULL; // object with which ray intersects first
    int index;
    for(int i = 0; i < objects.size(); i++){
        if(in_out == IN || exc != objects[i]){
            int num_intersects = objects[i] -> intersect_t(*r, intersects, normals);
            //printf("t: %f\n", t);
            if(num_intersects > 0 && intersects[0] < t && intersects[0] > 0.0){
                t = intersects[0];
                normal = normals[0];
                o = objects[i];
                index = i;
            }
        }
    }
    if(in_out == IN){
        printf("index: %d\n", index);
        printf("exc: %p, o: %p\n", exc, o);
        assert(exc == o);
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

GLdouble phong(point4d *p, vector4d *n, material *m, vector4d *in, vector4d *light, color *light_rgb, int channel){
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
    return m -> a.rgb[channel] * ambi.rgb[channel] +
        m -> d.rgb[channel] * (max_d(2, 0.0, normal * *light) * light_rgb -> rgb[channel]);
        //m -> s * (pow((max_d(2, 0.0, normal * half)) , m -> h) * *light_rgb);
}



GLdouble cast(ray *r, int channel){
    point4d intersect;
    vector4d normal;
    material m;
    firstHit(NULL, r, &intersect, &normal, &m, OUT);
    if(intersect.w == 0.0)
        return 0;
    GLdouble result = phong(&intersect, &normal, &m, r -> dir, &lsource, &source, channel);
    normal = ~normal;
    //color result = {normal.x, normal.y, normal.z};
    //assert_color(result);
    return result;
}
/*
color trace(ray* r, Object3D *prev, int depth) {
    ray flec, frac;
    color spec, refr, dull;
    color intensity, result;
    if(depth >= MAX_DEPTH){
        return back; 
    }
    point4d intersect = { 0, 0, 0, 0 };
    vector4d n;
    material m;
    Object3D *obj = firstHit(prev, r, &intersect, &n, &m);
    //n = ~n;
    if(intersect.w == 0){ // no intersect
        if((~(*(r->dir)) * ~lsource) == 1.0){
            result = source;
        }else
            result = back;
    } else {

        //printf("normal: ");
        //print_vector(n);
        // compute reflection
        if (m.s.r + m.s.g + m.s.b > 0) {
            vector4d norm = (*r->dir) | n;
            flec.copy(intersect, norm);
            color incoming = trace(&flec, obj, depth + 1);
            intensity = m.s * incoming;
            //double shine = pow(~(*r -> dir) * ~eye, m.h);
            //printf("shine: %f\n", shine);
            //intensity = (m.s * incoming + (color){shine, shine, shine}) * .5;
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
        point4d shadow_sect;
        firstHit(obj, &shadow, &shadow_sect, (vector4d*)nullptr, (material*)nullptr);
        if (shadow_sect.w == 0) {
            dull = phong(&intersect, &n, &m, r -> dir, &lsource, &source); 
            // no shadow (no intersection)
        } else {
            dull = phong(&intersect, &n, &m, r -> dir, &lsource, &back); 
        }
        result = (intensity + dull);
    }
    normalize(&result);
    return result;
}
*/

GLdouble trace(ray* r, Object3D *prev, int depth, int channel, int in_out) {
    ray flec, frac, shadow;
    GLdouble spec, refr, dull;
    GLdouble intensity, result;
    if(depth >= MAX_DEPTH){
        //printf("max depth\n");
        return back.rgb[channel]; 
    }
    point4d intersect = { 0, 0, 0, 0 };
    vector4d n;
    material m;
    Object3D *obj = firstHit(prev, r, &intersect, &n, &m, in_out);
    //n = ~n;
    if(intersect.w == 0){ // no intersect
        //printf("no intersect\n");
        if(((~(*r->dir)) * (~lsource)) == 1.0){
            result = source.rgb[channel];
        }else
            result = back.rgb[channel];
    } else {
        double refl, trans;
        double n1, n2;
        switch(in_out){
            case IN:
                n1 = m.b + m.c / w_length[channel];
                n2 = 1;
                break;
            case OUT:
                n1 = 1;
                n2 = m.b + m.c / w_length[channel];
                break;
        }
        //printf("n1: %f, n2: %f\n", n1, n2);
        double r0 = (n1 - n2) * (n1 - n2) / ((n1 + n2) * (n1 + n2));
        //printf("r0: %f\n", r0);
        vector4d norm = (*r->dir) | n; // reflection vector
        double cosi = (~norm) * (~n);
        double sini = sqrt(1 - cosi * cosi);
        double sint = sini * n1 / n2;
        double cost = sqrt(1 - sint * sint);
        //printf("cosi: %f, sini: %f, cost: %f, sint: %f\n", cosi, sini, cost, sint);
        bool tir = sini >= n2/n1; // total internal reflection
        if(n1 <= n2){
            refl = r0 + (1 - r0) * pow(1 - cosi, 5);
        }else if(n1 > n2 && !tir){
            refl = r0 + (1 - r0) * pow(1 - cost, 5);
        }else{
            refl = 1.0; 
        }
        if (m.s.rgb[channel] > 0) {
            //if(refl == 1.0)
            //    printf("reflectance: %f\n", refl);
            flec.copy(intersect, norm);
            spec = refl * m.s.rgb[channel] * trace(&flec, obj, depth + 1, channel, in_out);
            //double shine = pow(~(*r -> dir) * ~eye, m.h);
            //printf("shine: %f\n", shine);
            //intensity = (m.s * incoming + (color){shine, shine, shine}) * .5;
        } else {
            spec = 0;
        }

        // compute transmittance: 1 - reflectance
        trans = 1 - refl;
        double n0 = n1 / n2;
        vector4d t = n0 * (*r -> dir) + ((n0 * cosi - sqrt(1 - n0 * n0 * (1 - cosi *  cosi))) * n);
        //assert(t * *r -> dir >= 0.0);
        if (trans > 0 && !tir) {
            frac.copy(intersect, t);
            refr = trans * (1 - m.s.rgb[channel]) * trace(&frac, obj, depth + 1, channel, in_out ^ 1); 
            //printf("transmittance: %f\n", trans);
            if(obj == objects[7])
                printf("refracted: %f\n", refr);
            //refr = 0;
        } else {
            refr = 0;
        }
        // compute shadow
        shadow.copy(intersect, lsource);
        point4d shadow_sect;
        firstHit(obj, &shadow, &shadow_sect, (vector4d*)nullptr, (material*)nullptr, in_out);
        if (shadow_sect.w == 0) {
            dull = phong(&intersect, &n, &m, r -> dir, &lsource, &source, channel); 
            // no shadow (no intersection)
        } else {
            //printf("blocked\n");
            dull = phong(&intersect, &n, &m, r -> dir, &lsource, &back, channel); 
        }
        result = spec + refr + dull;
    }
    return result;
}
