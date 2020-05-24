.PHONY: sesc.opt sesc.debug

sesc.opt:
	scons -j 4 build/SMP_BOOKSIM/sesc.opt

sesc.debug:
	scons -j 4 build/SMP_BOOKSIM/sesc.debug
