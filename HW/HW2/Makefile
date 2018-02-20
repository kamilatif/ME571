CC = mpicc

SRC = ../../src
INCLUDE_DIR = -I../../src -I/usr/include

objs = ${SRC}/demo_util.c

.PHONY: p2 p3e p3b
 
p3b: p3b.c ${objs}
	${CC} -o p3b  $(INCLUDE_DIR) p3b.c ${objs} -lm

p3e: p3e.c ${objs}
	${CC} -o p3e  $(INCLUDE_DIR) p3e.c ${objs} -lm

p2: p2.c ${objs}
	${CC} -o p2  $(INCLUDE_DIR) p2.c ${objs} -lm
