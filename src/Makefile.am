# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${CPPFLAGS} ${PROJ_INCLUDES} ${PROJ_INCLUDES}
AM_LDFLAGS = ${LDFLAGS} ${PROJ_LDFLAGS} ${LIBS} ${PROJ_LDFLAGS} -lm


TARGETS = libuwpkfcvm.a libuwpkfcvm.so uwpkfcvm_query build_kdtree query_nearest_kdtree

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

kdtree_util.o: kdtree_util.c
	$(CC) $(AM_CFLAGS) -o $@ -c $^

build_kdtree.o: build_kdtree.c
	$(CC) $(AM_CFLAGS) -o $@ -c $^

build_kdtree: build_kdtree.o kdtree_util.o
	$(CC) -o $@ $^ $(AM_LDFLAGS)

query_nearest_kdtree.o: query_nearest_kdtree.c
	$(CC) $(AM_CFLAGS) -o $@ -c $^

query_nearest_kdtree: query_nearest_kdtree.o kdtree_util.o
	$(CC) -o $@ $^ $(AM_LDFLAGS)

clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

