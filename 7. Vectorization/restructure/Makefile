CC = cc
CFLAGS = -qopt-report=5 -restrict

OBJECTS = vectorisation.o

.SUFFIXES: .o .c

.c.o:
	$(CC) -c $(CFLAGS) -o "$@" "$<"

all: vect

vect: $(OBJECTS)
	$(CC) $(CFLAGS) -o vect $(OBJECTS)

clean: 
	rm -f *.o vect *.optrpt
