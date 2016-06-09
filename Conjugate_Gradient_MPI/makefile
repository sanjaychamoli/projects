CC = mpic++

SRC = cg.cpp
OBJ = $(SRC:.cpp=.o)

DEBUGSRC = Debug.cpp
DEBUGOBJ = $(DEBUGSRC:.cpp=.o)

ARRAYSRC = Array.cpp
ARRAYOBJ = $(ARRAYSRC:.cpp=.o)

DEBUG = -g
CFLAGS = $(DEBUG) -Wall -c -std=c++11 -Werror -Wshadow -Wextra
FLAGS = $(DEBUG) -Wall -std=c++11 -Werror -Wshadow -Wextra

PROGRAM = cg

All: $(PROGRAM)

$(PROGRAM) : $(DEBUGOBJ) $(SOLVEROBJ) $(OBJ)
	$(CC) $(FLAGS) $(DEBUGOBJ) $(SOLVEROBJ) $(OBJ) -o $(PROGRAM)

$(OBJ) : $(SOLVEROBJ) $(SRC)
	$(CC) $(CFLAGS) $(SRC)

$(DEBUGOBJ) : $(DEBUGSRC)
	$(CC) $(CFLAGS) $(DEBUGSRC)

$(SOLVEROBJ) : $(DEBUGOBJ) $(SOLVERSRC)
	$(CC) $(CFLAGS) $(SOLVERSRC)

clean:
	\rm *.o *~ $(PROGRAM) solution.txt exact.txt
