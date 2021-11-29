CC=gcc -g -Wall -lpthread -lm -Ilib
OBJECTS=raytracer_serial raytracer_threaded

all: $(OBJECTS)

$(OBJECTS): %: src/%.c
	$(CC) -o bin/$@ $<

clean::
	/bin/rm $(OBJECTS:%=./bin/%)
