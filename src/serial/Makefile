include ../../Makefile.in

# override for debugging
#COMPOPT=-g -O0 -DDEBUG

QSLE: QSLE.o qsle_main.o
	$(CPP) $(COMPOPT) $(II) $(LL) QSLE.o qsle_main.o -o ../../qsle_serial $(LIBS)

QSLE.o: QSLE.cpp ./include/QSLE.h
	$(CPP) $(COMPOPT) $(II) $(LL) -c QSLE.cpp $(LIBS)
qsle_main.o: qsle_main.cpp
	$(CPP) $(COMPOPT) $(II) $(LL) -c qsle_main.cpp $(LIBS)


clean:
	rm -rf ./*.o 
