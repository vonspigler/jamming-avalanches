SRC = SRC/
INC = SRC/INC/
CMP = COMP/
OPT_BEG = -std=c++11 -Wall -Wextra
OPT_END = -O3 #-lgsl -lgslcblas -lm
CXX = g++

all: MAIN.bin

MAIN.bin:\
  COMP/library.o COMP/system.o COMP/system_implementation.o COMP/initialize.o COMP/simulation_tools.o COMP/MAIN.o\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) COMP/MAIN.o COMP/simulation_tools.o COMP/initialize.o COMP/system_implementation.o COMP/system.o COMP/library.o -o $@ $(OPT_END)

COMP/MAIN.o:\
  $(SRC)MAIN.cpp\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) $(SRC)MAIN.cpp -c -o $@ $(OPT_END)

COMP/simulation_tools.o:\
  $(INC)simulation_tools.cpp\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) $(INC)simulation_tools.cpp -c -o $@ $(OPT_END)

COMP/initialize.o:\
  $(INC)initialize.cpp\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) $(INC)initialize.cpp -c -o $@ $(OPT_END)

COMP/system_implementation.o:\
  $(INC)system_implementation.cpp\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) $(INC)system_implementation.cpp -c -o $@ $(OPT_END)

COMP/system.o:\
  $(INC)system.cpp\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) $(INC)system.cpp -c -o $@ $(OPT_END)

COMP/library.o:\
  $(INC)library.cpp\
    $(INC)simulation_tools.h $(INC)initialize.h $(INC)system_implementation.h $(INC)system.h $(INC)library.h;\
  $(CXX) $(OPT_BEG) $(INC)library.cpp -c -o $@ $(OPT_END)


.PHONY: clean
clean:;\
  rm -f *~ MAIN.bin *.o fit.log;\
  rm -f COMP/*
