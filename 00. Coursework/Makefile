SRC=MD.c control.c
OBJ=$(SRC:.c=.o)
CC=cc
CFLAGS= -fast

all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -opt-report3 -lm

output.dat: MD input.dat
	./MD

clean:
	rm -f MD $(OBJ) output.dat* bench_c.e*

$(OBJ) : coord.h Makefile
