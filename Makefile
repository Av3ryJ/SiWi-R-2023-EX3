CXX = mpicxx
CXXFLAGS = -std=c++11 -Wall

INCLUDES =
LDFLAGS =
LIBS =

TARGET =
OBJS = $(TARGET).o

all: $(TARGET)

$(TARGET): $(TARGET).cpp Makefile
	$(CXX) $(CXXFLAGS) $< -o $(TARGET)


clean:
	@$(RM) -rf  $(TARGET)

