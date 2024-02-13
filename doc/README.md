# manysketch, fastgather, fastmultigather, multisearch, and manysearch - an introduction

This repository implements five sourmash plugins, `manysketch`, `fastgather`, `fastmultigather`, `multisearch`, and `manysearch`. These plugins make use of multithreading in Rust to provide very fast implementations of `sketch`, `search`, and `gather`. With large databases, these commands can be hundreds to thousands of times faster, and 10-50x lower memory, than sourmash.

The main *drawback* to these plugin commands is that their inputs and outputs are not as rich as the native sourmash commands. This may mean that your input files need to be prepared differently. The output may currently be most useful as a prefilter in conjunction with regular sourmash commands - see the instructions below for using `fastgather` to create picklists for sourmash.

## Input file formats


All four search/gather commands accept zip files, manifest files, or _text files containing lists of signature files_ ("fromfiles") for the search database. `multisearch`, `manysearch` and `fastmultigather` also use either zips, manifests, or "fromfiles" for queries, too. All commands now accept single signature files as well, though this is only useful for single-query input.

`manysketch` takes as input a CSV file with columns `name,genome_filename,protein_filename`. If you don't have `protein_filename` entries, be sure to include the trailing comma so the CSV reader can process the file correctly.

### Using zip files or manifest files


Manifest files are csv files with all information about sourmash signature parameters. Having a manifest allows us to select sketches relevant to the search (e.g. by k-mer size, scaled factor, etc) and perform checks without loading the sketches themselves into memory. We then only load the actual sketches (and optionally, downsample to a lower scaled value) when we're ready to use them.

If you have a `sourmash` zip file of signatures, it already contains a manifest that we can use internally.

If you'd like to generate a standalone `manifest` file from your signatures, you can do it like so:

```
sourmash sig manifest <sigs> -o sigs.manifest.csv
```
> Here, `sigs` can be any type of sourmash input, including a signature file or `pathlist`

### Using "fromfiles"

To prepare a **signature** fromfile from a database, first you need to split the database into individual files:
```
mkdir gtdb-reps-rs214-k21/
cd gtdb-reps-rs214-k21/
sourmash sig split -k 21 /group/ctbrowngrp/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k21.zip -E .sig.gz
cd ..
```

and then build a "fromfile":
```
find gtdb-reps-rs214-k21/ -name "*.sig.gz" -type f > list.gtdb-reps-rs214-k21.txt
```

When using these files for search, we have no a priori information about the parameters used for each sketch, so we load all signatures into memory at the start in order to generate a manifest. To avoid memory issues, the signatures are not kept in memory, but instead re-loaded as described below for each command (see: Notes on concurrency and efficiency). This makes using `pathlists` less efficient than `zip` files or `manifests`.

## Running the commands

### Running `manysketch`

The `manysketch` command sketches one or more FASTA/FASTQ files into a zipped sourmash signature collection (`zip`). `manysketch` uses one thread per input file, so it can (very) efficiently sketch many files at once; and, because sequence file parsing is entirely implemented in Rust, it is much faster than `sourmash sketch` for large FASTQ files.

To run `manysketch`, you need to build a text file list of FASTA/FASTQ files, with one on each line (`manysketch.csv`, below).  A simple way to do this for a directory is this command snippet:
```
echo name,genome_filename,protein_filename > manysketch.csv
for i in *.fa.gz
do
echo $i,$i,
done >> manysketch.csv
```


You can then run:

```
sourmash scripts manysketch manysketch.csv -o fa.zip
```
The output will be written to `fa.zip`

You can check if all signatures were written properly with
```
sourmash sig summarize fa.zip
```

To modify sketching parameters, use `--param-str` or `-p` and provide valid param string(s)
```
sourmash scripts manysketch fa.csv -o fa.zip -p k=21,k=31,k=51,scaled=1000,abund -p protein,k=10,scaled=200
```


### Running `multisearch`

The `multisearch` command compares one or more query genomes, and one or more subject genomes. It differs from `manysearch` by loading all genomes into memory.

`multisearch` takes two input collections (zip or "fromfiles"), and outputs a CSV:
```
sourmash scripts multisearch query-list.txt podar-ref-list.txt -o results.csv
```

To run it, you need to provide two collections of signature files. If you create a fromfile as above with GTDB reps, you can generate a query fromfile like so:

```
head -10 list.gtdb-reps-rs214-k21.txt > list.query.txt
```
and then run `multisearch` like so:

```
sourmash scripts multisearch list.query.txt list.gtdb-rs214-k21.txt  -o query.x.gtdb-reps.csv -k 21 --cores 4
```

The results file here, `query.x.gtdb-reps.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

### Running `fastgather`

The `fastgather` command is a much faster version of `sourmash gather`.

`fastgather` takes a single query metagenome (in any file format) and an input collection (zip or "fromfile") as database, and outputs a CSV:
```
sourmash scripts fastgather query.sig.gz podar-ref-list.txt -o results.csv --cores 4
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

Here the picklist should be used on a sourmash collection that contains a manifest - this will prevent sourmash from loading any sketches other than the ones in the fastgather CSV file. We recommend using zip file databases - manifests are produced automatically when `-o filename.zip` is used with `sketch dna`, and they also be prepared with `sourmash sig cat`. (If you are using a GTDB database, as above, then you already have a manifest!)

#### Example of picklist usage

A complete example Snakefile implementing the above workflow is available [in the 2023-swine-usda](https://github.com/ctb/2023-swine-usda/blob/main/Snakefile) repository. Note, it is slightly out of date at the moment!

### Running `fastmultigather`

`fastmultigather` takes a collection of query metagenomes and a collection of sketches as a database, and outputs many CSVs:
```
sourmash scripts fastmultigather query-list.txt podar-ref-lists.txt --cores 4
```

The main advantage that `fastmultigather` has over running `fastgather` on multiple queries is that you only load the database files once, which can be a significant time savings for large databases!

#### Output files for `fastmultigather`

`fastmultigather` will output two CSV files for each query, a `prefetch` file containing all overlapping matches between that query and the database, and a `gather` file containing the minimum metagenome cover for that query in the database.

The prefetch CSV will be named `{signame}.prefetch.csv`, and the gather CSV will be named `{signame}.gather.csv`.  Here, `{signame}` is the name of your sourmash signature.

**Warning:** At the moment, if two different queries have the same `{signame}`, the CSVs for one of the queries will be overwritten by the other query. The behavior here is undefined in practice, because of multithreading: we don't know what queries will be executed when or files will be written first.

### Running `manysearch`

The `manysearch` command compares one or more collections of query sketches, and one or more collections of subject sketches. It is the core command we use for searching petabase-scale databases of metagenomes for contained genomes.

`manysearch` takes two collections as input, and outputs a CSV:
```
sourmash scripts manysearch query-list.txt podar-ref-list.txt -o results.csv
```

To run it, you need to provide two "fromfiles" containing lists of paths to signature files (`.sig` or `.sig.gz`). If you create a fromfile as above with GTDB reps, you can generate a query fromfile like so:

```
head -10 list.gtdb-reps-rs214-k21.txt > list.query.txt
```
and then run `manysearch` like so:

```
sourmash scripts manysearch list.query.txt list.gtdb-rs214-k21.txt  -o query.x.gtdb-reps.csv -k 21 --cores 4
```

The results file here, `query.x.gtdb-reps.csv`, will have 8 columns: `query` and `query_md5`, `match` and `match_md5`, and `containment`, `jaccard`, `max_containment`, and `intersect_hashes`.

## Notes on concurrency and efficiency

Each command does things slightly differently, with implications for CPU and disk load. You can measure threading efficiency with `/usr/bin/time -v` on Linux systems, and disk load by number of complaints received when running.

(The below info is for fromfile lists. If you are using mastiff indexes, very different performance parameters apply. We will update here as we benchmark and improve!)

`manysketch` loads one sequence file from disk per thread and sketches it using all signature params simultaneously.

`manysearch` loads all the queries at the beginning, and then loads one database sketch from disk per thread. The compute-per-database-sketch is dominated by I/O. So your number of threads should be chosen with care for disk load. We typically limit it to `-c 32` for shared disks.

`multisearch` loads all the queries and database sketches once, at the beginning, and then uses multithreading to search across all matching sequences. For large databases it is extremely efficient at using all available cores. So 128 threads or more should work fine!

Like `multisearch`, `fastgather` loads everything at the beginning, and then uses multithreading to search across all matching sequences. For large databases it is extremely efficient at using all available cores. So 128 threads or more should work fine!

`fastmultigather` loads the entire database once, and then loads one query from disk per thread. The compute-per-query can be significant, though, so multithreading efficiency here is less dependent on I/O and the disk is less likely to be saturated with many threads. We suggest limiting threads to between 32 and 64 to decrease shared disk load.

## Appendix 1 - `index` to create a low-memory index

The command `sourmash scripts index` makes an on-disk inverted index
for low memory fast search. Indexing takes a while, but then search
takes fewer resources.

Currently only fastmultigather and manysearch can use this kind of index.
