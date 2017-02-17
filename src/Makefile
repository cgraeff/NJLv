TARGET = eos
LIB_DIR = -L./libdatafun/
LIBS = -lm -lgsl -lgslcblas -ldatafun
CC = gcc
CFLAGS = -Wall -Wno-unknown-pragmas

# Add path for gsl from Homebrew in OSX
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	INCLUDE = -I/usr/local/include
	LIB_DIR += -L/usr/local/lib
endif

.PHONY: default all clean

all: default

debug: CFLAGS += -g
debug: LIBFLAGS = debug
debug: exe

default: CFLAGS += -O3
default: exe

exe: datafun $(TARGET)

datafun:
	@cd libdatafun; make $(LIBFLAGS)

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	@$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	@$(CC) $(OBJECTS) -Wall $(LIB_DIR) $(LIBS) -o $@
	@mv $(TARGET) ..

clean:
	@-rm -f *.o
