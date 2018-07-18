##########################################################################
#
#  ===============================================
#   Schin's simple makefile  Ver. 990621
#  ===============================================
#
#   (Note)
#
#   (History)
#    990615 USERINC defined
#    990611 base1 introduced by S.D.
#    Jan-13-1997  R.Tanaka/SPring-8    created.
#
##########################################################################
#USERLIBS     =  -lcl -lM 
USERINC      = -I$(base1)/clib/head

DEFDIR       = $(base1)/cod2

OPTION       = -DDEBUG
#OPTION       = 

CFLAGS       = -e $(USERINC)
#CFLAGS       = -Ae $(USERINC)
######  -Ae option doesn't work on dec machines ######
######  -Ae option doesn't work on mac machines ######

LDFLAGS      =  $(USERLIBS) -lm

CC           =	cc  $(OPTION)

LINKER       =	cc  $(OPTION)
               
PROGRAM      = main

EXEIMAGE     = $(PROGRAM).exe

SOURCE       = $(PROGRAM).c

OBJS         = $(PROGRAM).o

#LOBJS        = $(base1)/clib/time/days.o\
#               $(base1)/clib/time/texprss.o
LOBJS        = 

all:         $(EXEIMAGE)

$(EXEIMAGE):  $(OBJS) $(LOBJS)
	     $(LINKER)  $(OBJS) $(LOBJS) -o $(EXEIMAGE) $(LDFLAGS)

$(OBJS):   $(SOURCE)
	$(CC) -c $(CFLAGS) $(SOURCE) -o $(OBJS)


clean:;      @rm -f $(OBJS) $(EXEIMAGE)

