SRC = \
 OpenETran.c \
 OEEngine.c \
 OERead.c \
 Parser.c \
 ReadUtils.c \
 WritePlotFile.c \
 ChangeTimeStep.c \
 Components/ArrBez.c \
 Components/Arrester.c \
 Components/BezUtils.c \
 Components/Capacitor.c \
 Components/Customer.c \
 Components/Ground.c \
 Components/Inductor.c \
 Components/Insulator.c \
 Components/Line.c \
 Components/LPM.c \
 Components/Meter.c \
 Components/Monitor.c \
 Components/NewArr.c \
 Components/PipeGap.c \
 Components/Pole.c \
 Components/Resistor.c \
 Components/Source.c \
 Components/SteepFront.c \
 Components/Surge.c \
 Components/Transformer.c

OBJ=$(SRC:.c=.o)

ifndef windir
RM=rm
EXE=OpenETran
else
RM=del
EXE=OpenETran.exe
endif

CC=gcc
CFLAGS=-Wall -O3 -I/gsl/gsl-2.1
LDFLAGS=-L/gsl/gsl-2.1/.libs -L/gsl/gsl-2.1/cblas/.libs -static -lgsl -lgslcblas -lm

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : all
all : $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@

.PHONY : clean
clean:
	-$(RM) $(OBJ)
