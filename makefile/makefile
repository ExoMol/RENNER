goal:   j-renner.x

tarball:
	tar cf renner.tar makefile *.f
        
checkin:
	ci -l Makefile *.f

###FOR  = ifort
FOR = ifort 
FFLAGS = -static -ipo   -openmp -O3 


LAPACK = -mkl=parallel

LIB     =   $(LAPACK) 

%.o : %.f
	$(FOR) -c $(FFLAGS) $<



OBJ    = cgamma.o cintens.o clock.o contra.o isort.o mues.o out.o pack.o \
         pfunc.o rencatst.o rendiag.o renmai.o renopt.o renrea.o rentel.o rentgst.o \
         rentzer.o rndelks.o rninter.o rnmopar.o rnmors.o rnplmi.o rtrhofct.o \
         rtstrt.o rtvibram.o rtwavene.o rtyfunc.o rwdip.o spnorb.o stropt.o \
         subdip.o triale.o trvals.o

j-renner.x:	$(OBJ)
	$(FOR) -o renner.x $(FFLAGS) $(OBJ) $(LIB)
clean:
	rm $(OBJ)

cgamma.o:	cgamma.f;	$(FOR) -c $(FFLAGS) cgamma.f
cintens.o:	cintens.f;	$(FOR) -c $(FFLAGS) cintens.f
clock.o:	clock.f;	$(FOR) -c $(FFLAGS) clock.f
contra.o:	contra.f;	$(FOR) -c $(FFLAGS) contra.f
isort.o:	isort.f;	$(FOR) -c $(FFLAGS) isort.f
mues.o:		mues.f;		$(FOR) -c $(FFLAGS) mues.f
out.o:		out.f;		$(FOR) -c $(FFLAGS) out.f
pack.o:		pack.f;		$(FOR) -c $(FFLAGS) pack.f
rencatst.o:	rencatst.f;	$(FOR) -c $(FFLAGS) rencatst.f
rendiag.o:	rendiag.f;	$(FOR) -c $(FFLAGS) rendiag.f
renmai.o:	renmai.f;	$(FOR) -c $(FFLAGS) renmai.f
pfunc.o:	pfunc.f;	$(FOR) -c $(FFLAGS) pfunc.f
renopt.o:	renopt.f;	$(FOR) -c $(FFLAGS) renopt.f
renrea.o:	renrea.f;	$(FOR) -c $(FFLAGS) renrea.f
rentel.o:	rentel.f;	$(FOR) -c $(FFLAGS) rentel.f
rentgst.o:	rentgst.f;	$(FOR) -c $(FFLAGS) rentgst.f
rentzer.o:	rentzer.f;	$(FOR) -c $(FFLAGS) rentzer.f
rndelks.o:	rndelks.f;	$(FOR) -c $(FFLAGS) rndelks.f
rninter.o:	rninter.f;	$(FOR) -c $(FFLAGS) rninter.f
rnmopar.o:	rnmopar.f;	$(FOR) -c $(FFLAGS) rnmopar.f
rnmors.o:	rnmors.f;	$(FOR) -c $(FFLAGS) rnmors.f
rnplmi.o:	rnplmi.f;	$(FOR) -c $(FFLAGS) rnplmi.f
rtrhofct.o:	rtrhofct.f;	$(FOR) -c $(FFLAGS) rtrhofct.f
rtstrt.o:	rtstrt.f;	$(FOR) -c $(FFLAGS) rtstrt.f
rtvibram.o:	rtvibram.f;	$(FOR) -c $(FFLAGS) rtvibram.f
rtwavene.o:	rtwavene.f;	$(FOR) -c $(FFLAGS) rtwavene.f
rtyfunc.o:	rtyfunc.f;	$(FOR) -c $(FFLAGS) rtyfunc.f
rwdip.o:	rwdip.f;	$(FOR) -c $(FFLAGS) rwdip.f
spnorb.o:	spnorb.f;	$(FOR) -c $(FFLAGS) spnorb.f
stropt.o:	stropt.f;	$(FOR) -c $(FFLAGS) stropt.f
subdip.o:	subdip.f;	$(FOR) -c $(FFLAGS) subdip.f
triale.o:	triale.f;	$(FOR) -c $(FFLAGS) triale.f
trvals.o:	trvals.f;	$(FOR) -c $(FFLAGS) trvals.f
