CC      =  gcc
CFLAGS  =  -g -Wall -D MAIN  #-O2  
LDFLAGS = -lz

PROG = sxspar 

.SUFFIXS:.c .o

all:$(PROG)

sxspar: sxs_reader.o sxs_writer.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) 	

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	


