Data = subCh* amr*

move:
	mkdir temp
	cp -f flash4 temp/
	cp -f flash.par temp/
	cp -f ofp_flash.sh temp/
	cp -f init.dat temp/
	mv -f ${Data} temp/

copy:
	cp ~/FLASH4.3/subCh/flash4 .

clean:
	rm -f *~

distclean:
	make clean
	rm -f ${Data} core
