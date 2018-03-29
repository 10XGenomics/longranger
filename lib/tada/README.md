# TDA_rust
Low-memory De Bruijn graph construction & path compression libraries.


## Roadmap

* move FASTQ handling & sorting to separate crate
* move martian stages to separate crate
* move other extraneous stuff out


### Kmer structs
* try to design a single-element kmer class that is generic over the integer types (u8, u16, u32, u64, u128) 
  * use the num crate to get a u128 implementation
  * when there is support for integer type parameters, use that to get intermediate sized kmers.
* define a kmer trait with the common kmer operations that can be reused
* update DeBruijn methods to accept the kmer trait
* build an AnnotatedKmer<K, T> struct that lets you attach generic data to a kmer, and proxies it's kmer implementation through to a kmer
* other option is to use the newtype pattern somehow to let the user define new structs that contain a kmer & implement the kmer trait. (would need to implement Deref?)




### DBG classes
* think about stranded / non-stranded analysis: RNA-seq generally knows strand - how to generalize for both cases?
* improve naming & ergonomics: TempGraph, Edge, VEdge, etc
* convenience methods for creating DBG without worrying about Lmer sharding
* path compression: can you pass a reducer in to summarize the per-kmer annotation data & associate with graph?
* path compression: can you provide a predicate on pairs of colors that controls where the path compression can start and stop

* think about APIs required for other operations on graphs -- bubble popping, edge trimming, etc.
