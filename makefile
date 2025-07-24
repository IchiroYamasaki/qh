# -----------------------------
# | makefile |
# -----------------------------

### macros ###
FC = ifort
LIBS = -framework Accelerate -Wl,-ld_classic, -heap-arrays 10
#DEBUG = -C                                                  
#DEBUG = -check all -traceback -g  

### file ###
TARGET1 = qh.exe
TARGET2 = test.exe
TARGET3 = qh_real.exe
TARGET4 = real.exe
TARGET5 = bdg_tw.exe
TARGET6 = bdg_proj.exe
TARGET7 = bdg.exe
TARGET8 = qh_tw.exe
TARGET9 = qh_tw_proj.exe

### rules for GNU make ###
%.o : %.f90
	$(FC) $(OPT) $(DEBUG) -c $< 
%.exe :
	$(FC) $^ $(LIBS) $(DEBUG) -o $@

### dependencies ###
$(TARGET1) : bacs.o functions.o LAP_ev.o ranpack.o Hamiltonian.o qh.o
$(TARGET2) : bacs.o ranpack.o functions.o Hamiltonian.o  test.o
$(TARGET3) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o qh_real.o
$(TARGET4) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o real.o
$(TARGET5) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o bdg_tw.o
$(TARGET6) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o bdg_proj.o
$(TARGET7) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o bdg.o
$(TARGET8) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o qh_tw.o
$(TARGET9) : bacs.o LAP_ev.o functions.o ranpack.o Hamiltonian.o qh_tw_proj.o

# time ./out
qh: $(TARGET1)
	time ./$(TARGET1)
	 make clean

test: $(TARGET2)
	time ./$(TARGET2)
	make clean

qh_real: $(TARGET3)
	time ./$(TARGET3)
	make clean

real: $(TARGET4)
	time ./$(TARGET4)
	make clean

bdg_tw: $(TARGET5)
	time ./$(TARGET5)
	make clean

bdg_proj: $(TARGET6)
	time ./$(TARGET6)
	make clean

bdg: $(TARGET7)
	time ./$(TARGET7)
	make clean

qh_tw: $(TARGET8)
	time ./$(TARGET8)
	make clean

qh_tw_proj: $(TARGET9)
	time ./$(TARGET9)
	make clean

# clean up
clean:
	rm -f *.o *.mod *~ *.exe
