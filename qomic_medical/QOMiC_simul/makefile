# makefile for creating vector and matrix library
 CC 			= gcc
CXX 			= g++ -std=c++11

LIBDIR			= ../Routines/
LIBDIR2			= ../Classes/

IFLAGS			= -I$(LIBDIR) -I/rds/general/user/aj1520/home/anaconda3/include/
LFLAGS			= -L$(LIBDIR) -L/rds/general/user/aj1520/home/anaconda3/lib/
LIBS			= -lgsl -lgslcblas
CXXFLAGS		= 
CPPFLAGS		= -W -Wall --pedantic -O3  $(IFLAGS)
#CPPFLAGS		= -W -Wall -pg -O0 -march=i686 $(IFLAGS)
#CPPFLAGS		= -W -Wall -O4 -g $(IFLAGS)

# Use the compile option -pg to time C++ code.

all:			QOMiC_simul .depend 

clean:		
		rm -f *.o
		rm -f .depend
		rm -f QOMiC_simul
		rm -f *~
		rm -f ../Routines/*.o
		rm -f ../Routines/*~
		rm -f ../Classes/*.o
		rm -f ../Classes/*~
		rm -f ../.depend
dep .depend:
	$(CXX) $(CPPFLAGS) -MM QOMiC_simul.cc >>.depend

infiles1		= QOMiC_simul.o  $(LIBDIR)dyn_name.o $(LIBDIR)matrix_handling.o  $(LIBDIR)rand.o $(LIBDIR2)Double_Matrices_cont.o $(LIBDIR2)Int_Matrices_cont.o $(LIBDIR)manip_dates.o $(LIBDIR)get_input.o $(LIBDIR)get_pbties.o

QOMiC_simul:	$(infiles1)
			$(CXX) $(CPPFLAGS) $(LFLAGS) -o QOMiC_simul $(infiles1) $(LIBS)

