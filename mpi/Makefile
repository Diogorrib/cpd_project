CC = gcc
CFLAGS = -O2 -fopenmp -Wall -Wextra -Werror #-g # uncomment to use vtune
LDFLAGS = -lm
SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)
DEPS = $(wildcard *.h)

TARGET = parsim

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
