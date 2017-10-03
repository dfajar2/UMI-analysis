CC = gcc
CFLAGS = -O2
LIBDIR = /sw/libs/bstring/cv/lib
INCDIR = /sw/libs/bstring/cv/include/bstring
	
fastq_qual_filter:	fastq_qual_filter.c 
	gcc -o fastq_qual_filter -I${INCDIR} -L${LIBDIR} -lbstring fastq_qual_filter.c 

clean:
	rm -f *.o fastq_qual_filter

