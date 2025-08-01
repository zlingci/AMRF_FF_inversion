FC=gfortran
optimize=-O
FFLAGS=$(optimize) -ffixed-line-length-none
SACHOME=/usr/local/sac
SACLIB= -L/$(SACHOME)/lib -lsac -lsacio -no-pie
CFLAGS=$(optimize)
CFLAGS=$(optimize) -DSAC_LIB

iterdeconfd.ye: iterdeconfd.ye.f recipes00.f bpfilter.f 
	$(FC) $(FFLAGS) -o iterdeconfd.ye bpfilter.f iterdeconfd.ye.f recipes00.f $(SACLIB) -lm

export:
	mkdir IterDecon
	cp iterdeconfd.ye.f IterDecon
	cp recipes00.f IterDecon
	cp bpfiter.f IterDecon
	cp Makefile_ye IterDecon
	tar cvf iterdecon.tar ./IterDecon
	compress iterdecon.tar
	rm -r IterDecon

clean:
	rm residual* d0 out.predict p1 numerator denominator 
