include ../Makefile.in

ALL: cubic_roots.so

cubic_roots.o: cubic_roots.c
	$(CC) -c $(CFLAGS) -fPIC cubic_roots.c -o cubic_roots.o

cubic_roots.so: cubic_roots.o
	$(LINK) $(LDFLAGS) cubic_roots.o -o cubic_roots.so

clean:
	rm -f *.o
	rm -f *.so
