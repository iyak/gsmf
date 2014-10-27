gsmf
====

Gibbs Sampling Motif Finder - simple implementation of randomized algorithm for motif finding [for study]

**Usage**

    $ make
    $ ./gsmf data/sample.fasta 4
    17	HOGE
    3	HOGE
    0	HOGE
    3	HOGE
    14	HOGE
    1	HOGE
    16	HOGE
    11	HOGE
    4	HOGE
    12	HOGE

first column shows positions of motif in each line.
detection sometimes fails as the nature of gibbs sampling.
