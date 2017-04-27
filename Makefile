SRC=MD.c control.c util.c
OBJ=$(SRC:.c=.o)
CC=cc
CFLAGS= -fast

all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -opt-report5 -lm

output.dat: MD input.dat
	./MD

clean:
	rm -f MD $(OBJ)

$(OBJ) : coord.h Makefile
