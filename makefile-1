CC = gcc
LD = gcc
CFLAGS = -O3 -Wall -Werror
LDFLAGS =-L/opt/X11/lib -lX11 -lm
INCLUDES=-I/opt/X11/include
RM = /bin/rm -f
OBJS = compare_gal_files.o file_operations.o graphics.o
EXECUTABLE = compare_gal_files

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE) $(LDFLAGS)

graphics.o: graphics.h graphics.c
	$(CC) $(CFLAGS) $(INCLUDES) -c graphics.c 

file_operations.o: file_operations.h file_operations.c
	$(CC) $(CFLAGS) -c file_operations.c

compare_gal_files.o: compare_gal_files.c graphics.h file_operations.h
	$(CC) $(CFLAGS) $(INCLUDES) -c compare_gal_files.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS)