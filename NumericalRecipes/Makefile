all: zrootsBinaryLens_wrapper.so

zrootsBinaryLens_wrapper.so: zrootsBinaryLens_wrapper.c complex.o laguer.o zroots.o nrutil.o
	gcc -lm -fPIC -Wl,-soname,rapper -shared -o $@ $? 

%.o: %.c
	gcc -lm -fPIC -c $^

.PHONY: clean
clean:
	$(RM) zrootsBinaryLens_wrapper.so complex.o laguer.o zroots.o nrutil.o

