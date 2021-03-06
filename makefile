CC = gcc

CFLAGS = -pg -g -ggdb -Wall -O3 -DNRAND -pthread
OBJS = color.o util.o colorant.o tabucol.o ant_fixed_k.o
BIN = colorant


all: $(OBJS)
	$(CC) -o $(BIN) $(CFLAGS) $(OBJS) -lm
	rm *.o

clean:
	rm $(OBJS) $(BIN) *~
