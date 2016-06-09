CXX = g++
CXXFLAGS = -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG

INCLUDES =
LDFLAGS =
LIBS =

# blas
INCLUDES += -I/usr/lib64/atlas/include/
LDFLAGS += -L/usr/lib64/atlas/
LIBS += -lcblas -latlas

# likwid
CXXFLAGS += -DUSE_LIKWID -pthread
INCLUDES += -I/usr/local/likwid-3.1.2/include/
LDFLAGS += -L/usr/local/likwid-3.1.2/lib/
LIBS += -llikwid

TARGET = matmult
OBJS = $(TARGET).o Array.o

all: $(TARGET)

$(TARGET): $(OBJS) Array.o Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) $(LIBS)

$(TARGET).o: $(TARGET).cpp multiplyfinal.h Array.h Timing.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(TARGET).cpp

Array.o: Array.cpp Array.h Makefile
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) Array.cpp

clean:
	@$(RM) -rf *.o $(TARGET)
