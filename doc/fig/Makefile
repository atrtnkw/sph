all:
	make -C shock_1d
	make -C shock_3d
	make -C strong
	make -C kh
	make -C pex
	make -C evrard

clean:
	rm -f *~

distclean:
	make clean
	make -C shock clean
	make -C strong clean
	make -C kh
	make -C pex clean
	make -C evrard clean
