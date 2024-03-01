# The rusty-neat project
NEAT (https://github.com/ncsa/neat) is essentially a simple program, since attempts to add some of the data packages from Python met with memory issues, time issues and more. As such, it required us to code NEAT in a way that would potentially translate well to a systems-level language such as C++ or Rust. Because of the desire to multi-thread NEAT, and the helpful way rust builds both thread parallelization and thread safety into their code, it is an appealing option to recreate NEAT and hopefully gain some speed and parallelization.

We will be taking an iterative approach, partly out of necessity as we learn the complexities of Rust. 

Check this branch for updates or submit a request to contribute code directly. There won't be any solid tasks yet until we're further along in development. 

How to use rusty-neat:

Download the executable in the release (current version 0.1.0).

```
./rusty-neat -h
```

displays help

```
./rusty-neat
```
rusty-neat on the default file (data/H1N1.fa). Only works from the repo level dir, run from anywhere else, you will need to specify the input reference file. The file is expected to be in FASTA format:

```angular2html
>Contig_name: other_info
AAAAAAAAAA
>Contig_name2: other_info2
GGGGGGGGGG
```

Use the help menu to see the available options and leave an issue if you find something bad happening. This data is not currently considered usable for anything requiring real rigor, but this is the first iteration toward a final product.
