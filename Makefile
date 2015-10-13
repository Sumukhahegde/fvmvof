NAME 	:= fvmvof
OBDIR	:= OBJ/
RUNDIR	:= RUN

OBJECTS =  $(OBDIR)particles.o   $(OBDIR)fvmvof.o 

HEADER 	:=	fvmvof.h\
		boundary.h

.DEFAULT:

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm 

CC	:=	gcc
FLAGS	:=	-g  $(INCPATH) -Wall -O3
#FLAGS	:=	-g -pg $(INCPATH) -Wall -fopenmp -O3
#LFLAGS	:=	$(LIBPATH) $(LIBS) -lm -lgsl -lblas
LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CC) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 
#	@cp iwrite $(RUNDIR)

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)particles.o: particles.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)particles.o particles.c

$(OBDIR)sph.o: sph.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)sph.o sph.c

copy:
	cp -f clean $(RUNDIR)/.	
clean: 		
	rm -f *~
	rm $(OBDIR)*.o
