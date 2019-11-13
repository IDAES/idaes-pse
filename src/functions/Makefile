include ../Makefile.in

ALL: functions.so

functions.o: functions.c
	$(CC) -c $(CFLAGS) functions.c -o functions.o

functions.so: functions.o
	$(CC) $(LDFLAGS) functions.o -o functions.so

clean:
	rm -f *.o
	rm -f *.so
