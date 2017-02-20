CC = gcc -std=c99 
LD = gcc -std=c99 
CFLAGS = -O3 -Wall -Werror
LDFLAGS =-L/opt/X11/lib -lX11 -lm 
INCLUDES=-I/opt/X11/include
RM = /bin/rm -f
OBJS = main.o file_operations.o graphics.o quad.o
EXECUTABLE = galsim

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE) $(LDFLAGS)

graphics.o: graphics.h graphics.c
	$(CC) $(CFLAGS) $(INCLUDES) -c graphics.c 

file_operations.o: file_operations.h file_operations.c
	$(CC) $(CFLAGS) -c file_operations.c

main.o: main.c graphics.h file_operations.h quad.h
	$(CC) $(CFLAGS) $(INCLUDES) -c main.c

quad.o: quad.h quad.c
	$(CC) $(CFLAGS) -c quad.c 

clean:
	$(RM) $(EXECUTABLE) $(OBJS)
