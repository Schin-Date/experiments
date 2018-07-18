##########################################################################
#
#  ===============================================
#   Standard Makefile for SR Software Development
#  ===============================================
#
#   (Note)
#    1. Do not change order of libraries and directories.
#    2. You can put your object and source files in OBJS or SRCS lines.
#
#   (History)
#    Jan-13-1997  R.Tanaka/SPring-8    created.
#
##########################################################################
# 4) Users
USERLIBS     =  -lcl -lM 

DEFDIR       = 

CFLAGS       = -Ae 

LDFLAGS      = $(USERLIBS) -lm

CC           =	cc 

LINKER       =	cc 
               
PROGRAM      = main

EXEIMAGE     = $(PROGRAM).exe

SOURCE       = $(PROGRAM).c

OBJS         = $(PROGRAM).o

LOBJS        = moon_pos.o\
               /users/schin/cod/celestial/read_nasa.o\
               /users/schin/cod2/lib/days.o\
               /users/schin/clib/stdio/inp_edt.o

all:         $(EXEIMAGE)

$(EXEIMAGE):  $(OBJS) $(LOBJS)
	     $(LINKER)  $(OBJS) $(LOBJS) -o $(EXEIMAGE) $(LDFLAGS)

$(OBJS):   $(SOURCE)
	$(CC) -c $(CFLAGS) $(SOURCE) -o $(OBJS)


clean:;      @rm -f $(OBJS) $(EXEIMAGE)

