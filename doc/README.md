# fastgather, fastmultigather, and manysearch - an introduction

This repository implements three sourmash plugins, `fastgather`, `fastmultigather`, and `manysearch`. These plugins make use of multithreading in Rust to provide very fast implementations of `search` and `gather`. With large databases, these commands can be hundreds to thousands of times faster, and 10-50x lower memory. 

The main *drawback* to these plugin commands is that their inputs and outputs are not as rich as the native sourmash commands. In particular, this means that input databases need to be prepared differently, and the output may be most useful as a prefilter in conjunction with regular sourmash commands.

## Preparing the database

All three commands use
_text files containing lists of signature files_, or "fromfiles", for the search database, and `manysearch` and `fastmultigather` use "fromfiles" for queries, too.

(Yes, this plugin will eventually be upgraded to support zip files; keep an eye on [sourmash#2230](https://github.com/sourmash-bio/sourmash/pull/2230).)

To prepare a fromfile from a database, first you need to split the database into individual files:
```
mkdir gtdb-reps-rs214-k21/
cd gtdb-reps-rs214-k21/
sourmash sig split -k 21 /group/ctbrowngrp/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k21.zip
cd ..
```

and then gzip the signatures and build a "fromfile":
```
find gtdb-reps-rs214-k21/ -name "*.sig" -exec gzip {} \;

find gtdb-reps-rs214-k21/ -type f > list.gtdb-reps-rs214-k21.txt
```
(Note that once [sourmash#2712](https://github.com/sourmash-bio/sourmash/pull/2712) is released, you will be able to do `sourmash sig split ... -E .sig.gz` to skip the additional find/gzip command.)

## Running the commands

### Running `manysearch`


The `manysearch` command finds overlaps between one or more query genomes, and one or more subject genomes.


`manysearch` takes two file lists as input, and outputs a CSV:
```
sourmash scripts manysearch query-list.txt podar-ref-list.txt -o results.csv
```

To run it, you need to provide two "fromfiles" containing lists of paths to signature files (`.sig` or `.sig.gz`). If you create a fromfile as above with GTDB reps, you can generate a query fromfile like so:

```
head -10 list.gtdb-reps-rs214-k21.txt > list.query.txt
```
and then run `manysearch` like so:

```
sourmash scripts manysearch list.query.txt list.gtdb-rs214-k21.txt  -o query.x.gtdb-reps.csv -k 21
```

The results file here, `query.x.gtdb-reps.csv`, will have five columns: `query` and `query_md5`, `match` and `match_md5`, and `containment.`

### Running `fastgather`

The `fastgather` command is a much faster version of sourmash gather.

`fastgather` takes a query metagenome and a file list as the database, and outputs a CSV:
```
sourmash scripts fastgather query.sig.gz podar-ref-list.txt -o results.csv
```

#### Using `fastgather` to create a picklist for `sourmash gather`

One handy use case for `fastgather` is to create a picklist that can be used by `sourmash gather`. This makes full use of the speed of `fastgather` while producing a complete set of `gather` outputs.

For example, if `list.gtdb-rs214-k21.txt` contains the paths to all GTDB RS214 genomes in `sig.gz` files, as above, then the following command will do a complete gather against GTDB:

```
sourmash scripts fastgather SRR606249.trim.sig.gz \
    list.gtdb-rs214-k21.txt -o SRR606249.fastgather.csv -k 21
```

This CSV file can then be used as a picklist for `sourmash gather` like so:

```
sourmash gather SRR606249.trim.sig.gz /group/ctbrowngrp/sourmash-db/gtdb-rs214/gtdb-rs214-k21.zip \
    --picklist SRR606249.fastgather.csv:match_name:ident \
    -o SRR606249.gather.csv
```

Here the picklist should be used on a sourmash collection that contains a manifest - this will prevent sourmash from loading any sketches other than the ones in the fastgather CSV file. We recommend zip files - zip file collections with manifests are produced automatically when `-o filename.zip` is used with `sketch dna`, and they also be prepared with `sourmash sig cat`. (If you are using a GTDB database, as above, then you already have a manifest!)

#### Example of picklist

A complete example Snakefile implementing the above workflow is available [in the 2023-swine-usda](https://github.com/ctb/2023-swine-usda/blob/main/Snakefile) repository.

### Running `fastmultigather`

`fastmultigather` takes a file list of query metagenomes and a file list for the database, and outputs many CSVs:
```
sourmash scripts fastmultigather query-list.txt podar-ref-lists.txt
```

The main advantage that `fastmultigather` has over `fastgather` is that you only load the database files once, which can be a significant time savings for large databases!

#### Output files for `fastmultigather`

`fastmultigather` will output two CSV files for each query, a `prefetch` file containing all overlapping matches between that query and the database, and a `gather` file containing the minimum metagenome cover for that query in the database.

The prefetch CSV will be named `{basename}.prefetch.csv`, and the gather CSV will be named `{basename}.gather.csv`.  Here, `{basename}` is the filename, stripped of its path.

**Warning:** At the moment, if two different queries have the same `{basename}`, the CSVs for one of the queries will be overwritten by the other query. The behavior here is undefined in practice, because of multithreading: we don't know what queries will be executed when or files will be written first.


