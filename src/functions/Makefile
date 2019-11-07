include ../Makefile.in

ifeq ($(ASL),)
    CFLAGS = -Wall -Wextra -O3 -fPIC
else
    CFLAGS = -Wall -Wextra -O3 -fPIC -I$(ASL)
endif

LDFLAGS = -shared -lm

ALL: functions.so

functions.o: functions.c
	$(CC) -c $(CFLAGS) functions.c -o functions.o

functions.so: functions.o
	$(CC) $(LDFLAGS) functions.o -o functions.so

clean:
	rm -f *.o
	rm -f *.so
