echidna
=======

:Author: Hyeshik Chang <hyeshik@snu.ac.kr>

``echidna`` is an on-the-fly stream parallelizer for FASTQ and FASTA
processing. It reads sequence data from standard input in either FASTQ
or FASTA format, and spreads them into several parallelized worker
subprocesses without creating temporary files.


Installation
------------

Follow the typical procedure of an autoconf-based package::

	./configure
	make
	make install

As the binary executable does not depend on any other file or installed
path, you can run it without installing.

Type ``./configure --help`` for more options.


Examples
--------

The following ``fastx_clipper`` command::

	zcat seq.fastq.gz | fastx_clipper -a CTGTAGGCACCATCAATC > seq-clipped.fq

Can be run parallel with ``echidna`` as following::

	zcat seq.fastq.gz | echidna fastx_clipper -a CTGTAGGCACCATCAATC > seq-clipped.fq

It runs as many ``fastx_clipper`` as cores in the machine.

Also, you can invoke pipelined commands like this::

	zcat seq.fastq.gz | fastq_quality_filter -q 30 -p 95 | \
		fastx_clipper -a CTGTAGGCACCATCAATC | fastq_to_fasta | gzip -c - > seq.fa.gz

With echidna for parallelization::

	zcat seq.fastq.gz | echidna -c "fastq_quality_filter -q 30 -p 95 | \
		fastx_clipper -a CTGTAGGCACCATCAATC | fastq_to_fasta" | gzip -c - > seq.fa.gz

The option ``-c`` is followed by a full shell command line as an argument.

You can change the number of processes to invoke with ``-p`` option::

	zcat seq.fastq.gz | echidna -p2 fastx_clipper -a CTGTAGGCACCATCAATC > seq-clipped.fq


Limitations
-----------

``Echidna`` distributes sequences into many subprocesses and collect them, it is
hard to keep the original order in the output.


Etymology
---------

The name is derived from "Echidna" in Greek mythology.
She is the mother of Ladon, Hydra, and Orthros, who all have many heads. :)

