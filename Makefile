###########################################################
# Project 1 Makefile

CC = g++
CFLAGS = -Wall -ggdb
INCLUDE = -I/lusr/X11/include -I/lusr/include
LIBDIR = -L/lusr/X11/lib -L/lusr/lib
# Libraries that use native graphics hardware --
# appropriate for Linux machines in Taylor basement
LIBS = -lglut -lGLU -lGL -lpthread -lm -std=c++11

###########################################################
# Options if compiling on Mac
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
CC = g++
CFLAGS = -Wall -g -D__MAC__ -std=c++11 -fdiagnostics-color=always
INCLUDE = 
LIBDIR = -L/usr/X11/lib
LIBS = -framework OpenGL -framework GLUT
endif

###########################################################
# Uncomment the following line if you are using Mesa
#LIBS = -lglut -lMesaGLU -lMesaGL -lm

raytrace: raytrace.cpp geometry.cpp light.cpp lowlevel.cpp raytrace.h lowlevel.h
	${CC} ${CFLAGS} ${INCLUDE} -o raytrace ${LIBDIR} raytrace.cpp geometry.cpp light.cpp lowlevel.cpp ${LIBS} 

clean:
	rm -f raytrace *.o core
