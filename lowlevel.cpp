/******************************************************************/
/*         Pixel operations                                       */
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

/* constants */

#define RED 0    /* offset to red byte */
#define GREEN 1    /* offset to green byte */
#define BLUE 2    /* offset to blue byte */

/* local data */

/* pointer to image */
GLubyte* canvas;

void initCanvas(int w, int h) {
  int i;

  /* get memory */
  canvas = (GLubyte*) malloc(w*h*3*sizeof(GLubyte));
  /* clear it */
  for (i=0; i<w*h*3; i++) {
    canvas[i]=0xFF;
  }
  width = w; 
  height = h;
}

void drawPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
  if ((x < 0) || (x >= width) || (y < 0) || (y >= height)) return;
  canvas[3*width*(y)+3*(x)+RED] = (char)(r*255);
  canvas[3*width*(y)+3*(x)+GREEN] = (char)(g*255);
  canvas[3*width*(y)+3*(x)+BLUE] = (char)(b*255);
}

/* draw the canvas array on the screen */
void flushCanvas() {
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glRasterPos3f(0.0,0.0,0.0);
  glDrawPixels(width,height,GL_RGB,GL_UNSIGNED_BYTE,canvas);
}

