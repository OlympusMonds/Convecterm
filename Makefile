TARGET = convecterm
LIBS = -lm -fopenmp
CC = gcc
CFLAGS = -Wall -Ofast -fopenmp -fopt-info-vec-missed #-DBENCH  # Add for benchmarking

.PHONY: default all clean

default: $(TARGET)
all: default

bench: clean all
	$(shell bench  './convecterm')

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
