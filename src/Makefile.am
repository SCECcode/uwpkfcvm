# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${PROJ_INCLUDES}
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS}


TARGETS = libuwpkfcvm.a libuwpkfcvm.so uwpkfcvm_query

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libuwpkfcvm.so ${prefix}/lib
	cp libuwpkfcvm.a ${prefix}/lib
	cp uwpkfcvm.h ${prefix}/include
	cp uwpkfcvm_query ${prefix}/bin

libuwpkfcvm.a: uwpkfcvm_static.o 
	$(AR) rcs $@ $^

libuwpkfcvm.so: uwpkfcvm.o
	$(CC) -shared $(AM_CFLAGS) -o libuwpkfcvm.so $^ $(AM_LDFLAGS)

uwpkfcvm.o: uwpkfcvm.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
uwpkfcvm_static.o: uwpkfcvm.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

uwpkfcvm_query.o: uwpkfcvm_query.c
	$(CC) $(AM_CFLAGS) -o $@ -c $^

uwpkfcvm_query : uwpkfcvm_query.o libuwpkfcvm.a
	$(CC) -o $@ $^ $(AM_LDFLAGS) -lm
	
clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

