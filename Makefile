CC=gcc -g -Wall -lpthread
CFILES=main.c
OFILES=$(CFILES:.c=.o)

jacobi:	$(OFILES)
	$(CC) -o main $(OFILES)

.c.o:
	$(CC) -c $<

clean::
	/bin/rm $(OFILES) jacobi
