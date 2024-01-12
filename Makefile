CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG

INCLUDES =
LDFLAGS =
LIBS =

TARGET = cg
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(TARGET).cpp Makefile
	$(CXX) $(CXXFLAGS) $< -o $(TARGET)


clean:
	rm -rf  $(TARGET)

