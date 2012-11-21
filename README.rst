echidna
=======

:Author: Hyeshik Chang <hyeshik@snu.ac.kr>

``echidna`` is an on-the-fly stream parallelizer for FASTQ and FASTA
processing. It reads sequence data from standard input in either FASTQ
or FASTA format, and spreads them into several parallelized worker
subprocesses without creating temporary files.
