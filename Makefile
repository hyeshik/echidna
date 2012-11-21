CFLAGS=	-ggdb -Wall -O3

all: echidna

echidna: echidna.c bufqueue.h
	${CC} -o $@ ${CFLAGS} echidna.c

longtest:
	date
	zcat longtest.fq.gz | ./echidna "fastx_clipper -a ATCAATTCGTATGCCGTCTTCTGCTTG" | pigz -p3 -c - > long.gz
	date

shorttest:
	zcat shorttest.gz|./echidna "fastx_clipper -a ATCAATTCGTATGCCGTCTTCTGCTTG" | pigz -p3 > shortresult.gz
	zcat shortresult.gz|wc -l
