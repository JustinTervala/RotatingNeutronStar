#/*************************************************************************
#*
#*                     MAKEFILE FOR QPO.C                                
#*                                                                         
#*************************************************************************/



#######################################################################
#                                                                     #
# 1. Specify C compiler and ANSI option:                              #
#                                                                     #      
####################################################################### 

#DEC ALPHA
#CC=cc -std1

#Linux
CC=g++ --std=c++11
COPTS=-Ofast -pipe -march=native -flto -funroll-loops -funsafe-math-optimizations -fno-omit-frame-pointer 
#/*************************************************************************
#*                     SPECIFY GRID SIZE
#*************************************************************************/

#STANDARD
SIZE=-DMDIV=65 -DSDIV=129

#HIGH
#SIZE=-DMDIV=101 -DSDIV=201

#VERY HIGH
#SIZE=-DMDIV=151 -DSDIV=301

#VERY VERY HIGH
#SIZE=-DMDIV=201 -DSDIV=401

#LOW
#SIZE=-DMDIV=51 -DSDIV=101

#VERY LOW
#SIZE=-DMDIV=41 -DSDIV=71

#/*************************************************************************
#*                     COMPILING FLAGS
#*************************************************************************/


# DEBUGGING OPTION
#MY_OWN =-g3

#/*************************************************************************
#*                    SOURCE AND OBJECT MACROS
#*************************************************************************/

OBJ=main.o RotatingNeutronStar.o EquationOfState.o GridTrig.o equil_util.o 

#/*************************************************************************
#*                    MAIN COMPILING INSTRUCTIONS
#*************************************************************************/

kepler: $(OBJ)
	$(CC) $(COPTS) $(MY_OWN) $(SIZE)  -o kepler $(OBJ) -lm -lrt


main.o: consts.h equil_util.h EquationOfState.hh RotatingNeutronStar.hh main.c
	$(CC)  $(COPTS) -c $(MY_OWN) $(CFLAGS) $(COPTFLAGS) $(SIZE)  main.c 

RotatingNeutronStar.o: RotatingNeutronStar.hh EquationOfState.hh consts.h equil_util.h GridTrig.hh
	$(CC) $(COPTS) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)  RotatingNeutronStar.cc

EquationOfState.o: EquationOfState.hh consts.h equil_util.h
	$(CC) $(COPTS) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)  EquationOfState.cc

equil_util.o:equil_util.h  consts.h equil_util.c
	$(CC)  $(COPTS) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)   equil_util.c

GridTrig.o: GridTrig.hh consts.h equil_util.h
	$(CC)  $(COPTS) -c $(MY_OWN) $(COPTFLAGS) $(SIZE)   GridTrig.cc
