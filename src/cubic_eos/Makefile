include ../Makefile.in

ifeq ($(ASL),)
    CFLAGS = -Wall -Wextra -O3
else
    CFLAGS = -Wall -Wextra -O3 -I$(ASL)
endif

LDFLAGS = -shared -lm

ALL: cubic_roots.so

cubic_roots.o: cubic_roots.c
	$(CC) -c $(CFLAGS) -fPIC cubic_roots.c -o cubic_roots.o

cubic_roots.so: cubic_roots.o
	$(CC) $(LDFLAGS) cubic_roots.o -o cubic_roots.so

clean:
	rm -f *.o
	rm -f *.so
