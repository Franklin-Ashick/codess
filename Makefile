CC=gcc
CFLAGS=-lm
ICC=icc

# Rule for making ci.exe
ci: cInsertion.c coordReader.c
	$(CC) cInsertion.c coordReader.c -o ci.exe $(CFLAGS)

# Rule for making fi.exe
fi: fInsertion.c coordReader.c
	$(CC) fInsertion.c coordReader.c -o fi.exe $(CFLAGS)

# Rule for making comp.exe with GNU compiler
comp: ompcInsertion.c coordReader.c
	$(CC) ompcInsertion.c coordReader.c -o comp.exe $(CFLAGS)

# Rule for making fomp.exe with GNU compiler
fomp: ompfInsertion.c coordReader.c
	$(CC) ompfInsertion.c coordReader.c -o fomp.exe $(CFLAGS)

# Rule for making icomp.exe with Intel compiler
icomp: ompcInsertion.c coordReader.c
	$(ICC) ompcInsertion.c coordReader.c -o icomp.exe

# Rule for making ifomp.exe with Intel compiler
ifomp: ompfInsertion.c coordReader.c
	$(ICC) ompfInsertion.c coordReader.c -o ifomp.exe

# 'make clean' to remove executables
clean:
	rm -f *.exe
