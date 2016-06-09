CC = g++
CFLAGS = -std=c++0x -Wall -Wextra -Wshadow -Werror -fopenmp -O3 -DNDEBUG 
RM = /bin/rm -f

SRC = main.cpp 
TARGET = rbgs
.PHONY : all clean depend
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $<

default: $(TARGET)
#	$(CC) -o $(TARGET) $(OBJS)

clean:
	$(RM) $(TARGET) solution.txt *.dat

depend:
	@makedepend -- $(CFLAGS) -- $(SRC) 
