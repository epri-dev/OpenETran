EXE=OpenETran.exe

SRC = \
 OpenETran.c \
 OEEngine.c \
 OERead.c \
 Parser.c \
 ReadUtils.c \
 WritePlotFile.c \
 ChangeTimeStep.c \
 components/Arrbez.c \
 components/Arrester.c \
 components/BezUtils.c \
 components/Capacitor.c \
 components/Customer.c \
 components/Ground.c \
 components/Inductor.c \
 components/Insulator.c \
 components/Line.c \
 components/LPM.c \
 components/Meter.c \
 components/Monitor.c \
 components/NewArr.c \
 components/PipeGap.c \
 components/Pole.c \
 components/Resistor.c \
 components/Source.c \
 components/SteepFront.c \
 components/Surge.c \
 components/Transformer.c

OBJ=$(SRC:.c=.o)

CC=gcc
CFLAGS=-Wall -O3 -I/gsl/gsl-2.1
LDFLAGS=-L/gsl/gsl-2.1/.libs -L/gsl/gsl-2.1/cblas/.libs -Bstatic -llibgsl -llibgslcblas
RM=rm

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : all
all : $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@

.PHONY : clean
clean:
	-$(RM) $(OBJ) core
