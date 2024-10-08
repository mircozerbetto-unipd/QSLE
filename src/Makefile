include ../Makefile.in

# override for debugging
#COMPOPT=-g -O0 -DDEBUG

QSLE: QSLE_parallel.o qsle_main.o
	$(CPP) $(COMPOPT) $(II) $(LL) QSLE_parallel.o qsle_main.o -o ../qsle $(LIBS)

QSLE_parallel.o: QSLE_parallel.cpp ./include/QSLE_parallel.h
	$(CPP) $(COMPOPT) $(II) $(LL) -c QSLE_parallel.cpp $(LIBS)
qsle_main.o: qsle_main.cpp
	$(CPP) $(COMPOPT) $(II) $(LL) -c qsle_main.cpp $(LIBS)


clean:
	rm -rf ./*.o
