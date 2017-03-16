compiler = gfortran -g -fbacktrace -ffpe-trap=zero,overflow,underflow

objs = retard.o fft.o plodia_venturia_model1.o

model: $(objs)
	$(compiler) -o model $(objs)

retard.o: retard.f
	$(compiler) -c retard.f
fft.o: fft.f
	$(compiler) -c fft.f
plodia_venturia_model1.o: plodia_venturia_model1.f
	$(compiler) -c plodia_venturia_model1.f

clear:
	rm -f *.o model

clean:
	rm -f results/*.data

run:
	if [ ! -d results ]; then mkdir results; fi
	./model
