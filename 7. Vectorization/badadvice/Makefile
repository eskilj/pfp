CC = cc
CFLAGS = -qopt-report=5

OBJECTS = external.o main.o

.SUFFIXES: .o .c

.c.o:
	$(CC) -c $(CFLAGS) -o "$@" "$<"

all: badprog

badprog: $(OBJECTS)
	$(CC) $(CFLAGS) -o badprog $(OBJECTS)

clean: 
	rm -f *.o badprog *.optrpt
