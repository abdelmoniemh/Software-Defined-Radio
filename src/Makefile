# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3 -g -std=c++20
# CFLAGS = -c
LDFLAGS =
INCLUDE = -I../include/
SOURCES = project.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp audioBlock.cpp PLL.cpp rdsUtil.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) -lpthread

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
