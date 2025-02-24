CXX = g++
CXXFLAGS = -O2 -fopenmp -Wall -Wextra -Werror
LDFLAGS = -lm

SRCS = $(wildcard *.cpp)  # Change *.c to *.cpp if your code is in C++
OBJS = $(SRCS:.cpp=.o)
DEPS = $(wildcard *.h)

TARGET = parsim

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
