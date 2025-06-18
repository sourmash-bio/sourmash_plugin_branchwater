# The branchwater plugin for sourmash

<p align="center"><img src="https://raw.githubusercontent.com/sourmash-bio/sourmash_plugin_branchwater/main/doc/_static/logo.png" height="256" /></p>

| command | functionality | docs |
| -------- | -------- | -------- |
| `manysketch` | Rapidly build sketches for many input files     | [link](#Running-manysketch)     |
| `singlesketch` | Sketch a single sample | [link](#Running-singlesketch)
| `fastgather` | Multithreaded `gather` of **one** metagenome against a database| [link](#Running-fastgather)
| `fastmultigather` | Multithreaded `gather` of **multiple** metagenomes against a database | [link](#Running-fastmultigather)
| `manysearch` | Multithreaded containment search for many query genomes in many large metagenomes | [link](#Running-manysearch)
| `multisearch` | Multithreaded comparison of multiple sketches, in memory | [link](#Running-multisearch-and-pairwise)
| `pairwise` | Multithreaded pairwise comparison of multiple sketches, in memory | [link](#Running-multisearch-and-pairwise)
| `cluster` | cluster sequences based on similarity data from `pairwise` or `multisearch` | [link](#Running-cluster)
| `index` | build a RocksDB inverted index for efficient containment queries | [link](#Running-index)

This repository implements multithreaded plugins for
[sourmash](https://sourmash.readthedocs.io/) that provide very fast
implementations of `sketch`, `search`, and `gather`. These commands
are typically hundreds to thousands of times faster, and 10-50x lower
memory, than the current sourmash code. For example, a `gather` of
SRR606249 with sourmash v4.8.6 against GTDB rs214 takes 40 minutes and
14 GB of RAM, while `fastgather` with 64 cores takes only 2 minutes
and 2 GB of RAM.

The main *drawback* to these plugin commands is that their inputs are
more restricted than the native sourmash commands, and in some cases
their outputs are in a different format. This means that your input
files may need to be prepared differently, and the output may need to
be processed differently.  The plugin commands are also a bit less
user friendly, because (for now) we're more focused on speed than
polish and user experience.

## Input file formats

sourmash supports a variety of different storage formats for sketches (see [sourmash docs](https://sourmash.readthedocs.io/en/latest/command-line.html#choosing-signature-output-formats)), and the branchwater plugin works with some (but not all) of them.

**We recommend using zip files or standalone manifest CSVs pointing to zip files whenever you need to provide multiple sketches.**

| command | command input | database format |
| -------- | -------- | -------- |
| `manysketch`     | CSV with input fasta/fastq paths (details below)    | _produces_ Zip database |
| `gather`     | Single metagenome in sig, zip, or pathlist     | Zip or pathlist |
| `fastmultigather` | Multiple metagenomes in sig, zip, or pathlist | Zip, pathlist, or rocksdb index |
| `manysearch` | Multiple genomes in sig, zip, or pathlist | Zip, pathlist, or rocksdb index |
| `multisearch` | Multiple sketches in sig, zip, or pathlist | Multiple sketches in sig, zip, or pathlist |
| `pairwise` | Multiple sketches in sig, zip, or pathlist | N/A |
| `cluster`| Output from `pairwise` or `multisearch`| N/A |
| `index` | Multiple sketches in sig, zip, or pathlist | N/A |

### Using zipfiles

When working with large collections of small sketches such as genomes, we suggest using zipfiles as produced by sourmash (e.g. using `sourmash sig cat` or `manysketch`). Zip files have a few nice features:

* sketches are compressed in zip files;
* zip files can contain many sketches, including incompatible types (e.g. multiple k-mer sizes);
* subsets of zip files can be efficiently selected and loaded;
* in particular, _single_ sketches can be loaded on demand, supporting lower memory requirements for certain kinds of searches.

For all these reasons, zip files are the most efficient and effective
basic storage type for sketches in sourmash, and the branchwater
plugin fully supports them!

You can create zipfiles with sourmash like so:
```
sourmash sig cat <list of sketches> -o sigs.zip
```

### Using manifests for input databases - why and when?

The branchwater plugin commands take a relatively restricted set of
inputs, compared to sourmash: they take a single file for query and a
single file for search.  These files must be a collection of sketches
(in .sig/.sig.gz format), a RocksDB database, a pathlist, or a
standalone manifest.

When you want to search subsets or deal with many large
files, we recommend _standalone manifests_ over constructing new
databases or using pathlists to point at subsets of files. Why?

There are three main reasons NOT to use zip files or pathlists in this case:

* first, metagenome sketches are often extremely large (100s of MBs to
  GBs), and it is not ideal to zip them all into a single
  zip file;
* second, both `manysearch` and `fastmultigather` take a single
  argument that specifies collections of metagenomes which need to be
  loaded on demand, because they cannot fit into memory;
* third, when searching subsets of large collections, it can be expensive
  to make a copy of the subset;

Manifests are a sourmash filetype that catalog sketch content and can
be used to point at many sketches and/or subset large collections of
sketches, and so they provide solutions to these problems. In
particular, manifests let you provide large collections of
sketches/collections of large sketches to `manysearch` and
`fastmultigather`, and also let you select a subset of a collection
without making a copy.

The branchwater plugin supports manifest CSVs.  These can be created
from lists of sketches by using `sourmash sig collect` or `sourmash
sig check`; for example,

```
sourmash sig check <database.sig.zip -F csv --picklist idents.csv:ident:ident -m subset.csv
```
will create a manifest CSV containing a subset of sketches [using picklists](https://sourmash.readthedocs.io/en/latest/command-line.html#using-picklists-to-subset-large-collections-of-signatures),

and

```
find -type f /path/to/sig/files/* > pathlist.txt
sourmash sig collect pathlist.txt -o summary-manifest.csv -F csv
```
will collect a list of all of the sketches under `/path/to/sig/files`
and make the list available through a combined manifest.

Note here that manifests are _much_ smaller than the files containing all
of the sketches!

Note also that manifests have many advantages over pathlists: in
particular, they contain metadata that enables fast loading of
specific sketches, and they support subsetting from large databases;
pathlists support neither.

### Using RocksDB inverted indexes

RocksDB indexes provide fast and low-memory lookups when searching
very large datasets, and are used for the branchwater petabase scale
search hosted at
[branchwater.sourmash.bio](https://branchwater.sourmash.bio).  As of
sourmash v4.9.0, RocksDB is fully supported by sourmash, as well.
See [the RocksDB HOWTO](https://sourmash.readthedocs.io/en/latest/howto-rocksdb.html) for more details.

Some commands - `fastmultigather` and `manysearch` - support using
these RocksDB-based inverted indexes for efficient search, and they
can be created by running `sourmash scripts index`. See
[the `index` documentation, below](#Running-index).

### Using "pathlists"

**Note: We no longer recommend using "pathlists". Use zip files or
  standalone manifests instead.**

You can make a pathlist by listing a collection of .sig.gz files like so:
```
find /path/to/directory/ -name "*.sig.gz" -type f > directory.txt
```

When using a pathlist for search, we load all signatures into memory
at the start in order to generate a manifest. To avoid memory issues,
the signatures are not kept in memory, but instead re-loaded as
described below for each command (see: Notes on performance, concurrency and
efficiency at the bottom). This makes using pathlists less efficient than `zip`
files.

## Running the commands

### Running `manysketch`

The `manysketch` command sketches a list of FASTA/FASTQ files into a
zipped sourmash signature collection (`zip`). `manysketch` uses one
thread per input file, so it can (very) efficiently sketch many files
at once; and, because sequence file parsing is entirely implemented in
Rust, it is much, _much_ faster than `sourmash sketch` for large FASTQ
files. However, it does not currently support translation,
i.e. protein signature generation from DNA FASTA.

`manysketch` can be used to sketch collections of genomes, pairs of R1/R2
files, and more flexible collections based on prefix - which approach is
used depends on the columns provided to `manysketch` in the input CSV file.

#### Specifying input FASTA

To run `manysketch`, you need to build a text file list of FASTA/FASTQ files (see `manysketch.csv` example, below).

The following formats are accepted:
- 3 columns: `name,genome_filename,protein_filename`
  >`genome_filename` entries are considered DNA FASTA, `protein_filename` entries are considered protein FASTA.
- 3 columns: `name,read1,read2`
  > All entries considered DNA FASTA, and both `read1` and `read2` files are used as input for a single sketch with name `name`.
- 4 columns: `name,input_moltype,prefix,exclude`
  > This filetype uses `glob` to find files that match `prefix` but do not match `exclude`. As such, `*` are ok in the `prefix` and `exclude` columns. Since we are dealing with "prefixes" here, we automatically search with `*` on the end of the `prefix` entry.

A simple way to build a manysketch input file for a directory is this command snippet:
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
See [the sourmash sketch docs](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-sketch-make-sourmash-signatures-from-sequence-data) for more information on param strings.

#### singleton sketching

`manysketch` also supports building independent sketches for each record in a FASTA file (`--singleton`).

You can run:

```
sourmash scripts manysketch manysketch.csv -o fa.zip --singleton
```
The output will be written to `fa.zip`

You can check if all signatures were written properly with
```
sourmash sig summarize fa.zip
```
The number of sketches per parameter combination should equal the total number of records in all input FASTA.
The `name` column will not be used. Instead, each sketch will be named from the FASTA record name.

#### Protein sketching: hp and dayhoff moltypes

`manysketch` supports all sourmash moltypes: `dna`, `skipm1n3`, `skipm2n3`, `protein`, `hp`, and `dayhoff`. See the `sourmash` [DNA](https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html#dna-encodings) and [protein](https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html#protein-encodings) encoding docs, as well as the [`sourmash` parameter documentation](https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html#default-parameters) for more information about what these "moltypes" mean and their default parameters.

`manysketch` does not translate DNA to protein, sorry. You'll need to do that ahead of time.

If you have a `proteins.csv` file which looks like:

```
name,genome_filename,protein_filename
protein1,,protein1.fa
protein2,,protein2.fa
```

You can run:

```
sourmash scripts manysketch proteins.csv -o proteins.zip -p dayhoff,k=16 -p protein,k=10 -p hp,k=42
```
The output will be written to `proteins.zip`

You can check if all signatures were written properly with

```
sourmash sig summarize proteins.zip
```

In this case, three sketches of `protein`, `dayhoff`, and `hp` moltypes were made for each file in `proteins.csv` and saved to `proteins.zip`.


## Running `singlesketch`

The `singlesketch` command generates a sketch for a single sample, with one or more input FASTA/FASTQ files.

### Basic Usage

```bash
sourmash scripts singlesketch input.fa -p k=21,scaled=1000,dna -o output.sig --name signature_name
```

### Using `stdin/stdout`

You can use `-` for `stdin` and output the result to `stdout`:
```bash
cat input.fa | sourmash scripts singlesketch - -o -
```

### Running `multisearch` and `pairwise`

The `multisearch` command compares one or more query genomes, and one or more subject genomes. It differs from `manysearch` because it loads everything into memory.

`multisearch` takes two input collections and outputs a CSV:
```
sourmash scripts multisearch query.sig.gz database.zip -o results.csv
```

The results file `results.csv`, will have 8 columns: `query` and
`query_md5`, `match` and `match_md5`, and `containment`, `jaccard`,
`max_containment`, and `intersect_hashes`.

The `pairwise` command does the same comparisons as `multisearch` but
takes only a single collection of sketches, for which it calculates
all the pairwise comparisons. Since the comparisons are symmetric, it
is approximately twice as fast as `multisearch`.

The `-t/--threshold` for `multisearch` and `pairwise` applies to the
containment of query-in-target and defaults to 0.01.  To report _any_
overlap between two sketches, set the threshold to 0.

### Running `fastgather`

The `fastgather` command is parallelized (and typically much faster)
version of `sourmash gather`.

`fastgather` takes a single query metagenome and a database, and outputs a CSV:
```
sourmash scripts fastgather query.sig.gz database.zip -o results.csv --cores 4
```

`fastgather` outputs the same columns as `sourmash gather`, with only a few exception
* `match_name` is output instead of `name`;
* `match_md5` is output instead of `md5`;
* `match_filename` is output instead of `filename`, and the value is different;
* `potential_false_negative` is not present in `fastgather` output;

### Running `fastmultigather`

`fastmultigather` takes a collection of query metagenomes and a collection of genome sketches as a database, and outputs a CSV file containing all the matches.
```
sourmash scripts fastmultigather queries.manifest.csv database.zip --cores 4 --save-matches -o results.csv
```

We suggest using standalone manifest CSVs wherever possible, especially if
the queries are large.

The main advantage that `fastmultigather` has over running
`fastgather` on multiple queries is that `fastmultigather` only needs
to load the database once for all queries, unlike with `fastgather`;
this can be a significant time savings for large databases.

#### Output files for `fastmultigather`

`fastmultigather` will output a gather file containing all results in
one file, specified with `-o/--output`. `fastmultigather` gather CSVs
provide the same columns as `fastgather`, above.

In addition, on a database of sketches (but not on RocksDB indexes)
`fastmultigather` will output a `prefetch` file containing all
overlapping matches between that query and the database. The prefetch
CSV will be named `{signame}.prefetch.csv`, where `{signame}` is the
name of your sourmash signature.

`--save-matches` is an optional flag that will save the matched hashes
for each query in a separate sourmash signature
`{signame}.matches.sig`. This can be useful for debugging or for
further analysis.

**Warning:** At the moment, if two different queries have the same
  `{signame}`, the output files for one query will be overwritten by
  the results from the other query. The behavior here is undefined in
  practice, because of multithreading: we don't know what queries will
  be executed when or files will be written first.

### Running `manysearch`

The `manysearch` command searches for query
sketches (genomes) in one or more collections of subject sketches (metagenomes). It is the
core command we use for searching petabase-scale databases of
metagenomes for contained genomes.

`manysearch` takes two collections as input, and outputs a CSV:
```
sourmash scripts manysearch queries.zip metagenomes.manifest.csv -o results.csv
```
We suggest using a manifest CSV for the metagenome collection.

The results file here, `query.x.gtdb-reps.csv`, will have the
following columns: `query`, `query_md5`, `match_name`, `match_md5`,
`containment`, `jaccard`, `max_containment`, `intersect_hashes`,
`query_containment_ani`.

If you run `manysearch` _without_ using a RocksDB database (that is,
against regular sketches), the results file will also have the
following columns: , `match_containment_ani`,
`average_containment_ani`, and `max_containment_ani`.

Finally, if using sketches that have abundance information, the
results file will also contain the following columns: `average_abund`,
`median_abund`, `std_abund`, `n_weighted_found`, and
`total_weighted_hashes`.

See
[the prefetch CSV output column documentation](https://sourmash.readthedocs.io/en/latest/classifying-signatures.html#appendix-e-prefetch-csv-output-columns)
for information on these various columns.

The `-t/--threshold` for `manysearch` applies to the containment of
query-in-database (e.g. genome in metagenome) and defaults to 0.01.
To report _any_ overlap between two sketches, set the threshold to 0.
(This will produce many, many results when searching a collection of
metagenomes!)

Using `-A/--output-all-comparisons` will ignore the threshold parameter
and output all comparisons done. Against a RocksDB database, only matches
with some overlap will be reported; with collections of sketches, all
pairs will be reported.

By default, `manysearch` will display the contents of the CSV file in a
human-readable format. This can be disabled with `-N/--no-pretty-print`
when executing large searches.

### Running `cluster`

The `cluster` command conducts graph-based clustering via the sequence
similarity measures in `pairwise` or `multisearch` outputs. It is a
new command and we are exploring its utility.

`cluster` takes the csv output of `pairwise` or `multisearch` input,
and outputs two CSVs:

1. `-o`, `--output` will contain the names of the clusters and the
   `ident` of each sequence included in the cluster
   (e.g. `Component_1, name1;name2`)

```
cluster,nodes
Component_1,name1;name2;name3
Component_2,name4
```

2. `--cluster-sizes` will contain information on cluster size, with a
   counts for the number of clusters of that size. For the two
   clusters above, the counts would look like this:

```
cluster_size,count
3,1
1,1
```

`cluster` takes a `--similarity_column` argument to specify which of
the similarity columns, with the following choices: `containment`,
`max_containment`, `jaccard`, `average_containment_ani`,
`maximum_containment_ani`. All values should be input as fractions
(e.g. 0.9 for 90%)

### Running `index`

The `index` subcommand creates a RocksDB inverted index that can be
used as an efficient database for `manysearch` (containment queries
into mixtures) and `fastmultigather` (mixture decomposition against a
database of genomes). RocksDB databases created by `index` in this
plugin are fully compatible with RocksDB databases created by
`sourmash index`, and vice versa.

RocksDB inverted indexes support fast, low-latency, and low-memory
queries, in exchange for a time-intensive indexing step and extra disk
space. They are also used for the
[branchwater.sourmash.bio](https://branchwater.sourmash.bio/)
real-time SRA metagenome query.

Note that RocksDB indexes do not support abundance estimation for
containment queries.

To run `index`, provide it with multiple sketches in sig, zip, or
pathlist format, and specify the desired output directory; we suggest
using the `.rocksdb` extension for RocksDB databases, e.g. `-o
gtdb-rs214-k31.rocksdb`.

By default, `index` will store a copy of the sketches
along with the inverted index.  This will substantially increase the
disk space required for large databases.  You can provide an optional
`--no-internal-storage` to `index` to store them externally, which
reduces the disk space needed for the index.  Read below for technical
details!

`index` can take any of the supported input types, but
unless you are using a zip file or a pathlist of JSON files, it may
need to load all the sketches into memory before indexing
them. Moreover, you can only use external storage with a zip file. We
are working on improving this; see
[issue #415](https://github.com/sourmash-bio/sourmash_plugin_branchwater/issues/415)
for details. A warning will be printed to stderr in situations where
the sketches are being loaded into memory.

#### Internal vs external storage of sketches in a RocksDB index

RocksDB indexes support containment queries (a la the
[branchwater application](https://github.com/sourmash-bio/branchwater)),
as well as `gather`-style mixture decomposition (see
[Irber et al., 2022](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2)).
For this plugin, the `manysearch` command supports a RocksDB index for
the database for containment queries, and `fastmultigather` can use a
RocksDB index for the database of genomes.

RocksDB indexes contain references to the sketches used to construct
the index. If `--internal-storage` is set (which is the default), a
copy of the sketches is stored within the RocksDB database directory;
if `--no-internal-storage` is provided, then the references point to
the original source sketches used to construct the database, wherever
they reside on your disk.

The sketches *are not used* by `manysearch`, but *are used* by
`fastmultigather`: you'll get an error if you
run `fastmultigather` against a RocksDB index where the sketches
cannot be loaded.

What this means is therefore a bit complicated, but boils down to
the following two approaches:

1. The safest thing to do is build a RocksDB index and use internal
storage (the default). This will consume more disk space but your
RocksDB database will always be usable for both `manysearch` and
`fastmultigather`, as well as the branchwater app.
2. If you want to avoid storing duplicate copies of your sketches,
then specify `--no-internal-storage` and provide a stable absolute
path to the source sketches.  This will again support both
`manysearch` and `fastmultigather`, as well as the branchwater app.
If the source sketches later become unavailable, `fastmultigather`
will stop working (although `manysearch` and the branchwater app
should be fine).

You should (for the moment) avoid specifying relative paths to the
sketches when running `index`.  Follow
[sourmash_branchwater_plugin#415](https://github.com/sourmash-bio/sourmash_plugin_branchwater/issues/415)
if better support for relative paths is of interest!

RocksDB databases created by `sourmash index` always use internal storage,
but RocksDB databases created with external storage using this plugin
can still be used by sourmash.

#### Links and more materials

Note that RocksDB indexes are implemented in the core
[sourmash Rust library](https://crates.io/crates/sourmash), and used
in downstream software packages (this plugin, and
[the branchwater application code](https://github.com/sourmash-bio/branchwater)).
The above documentation applies to sourmash core v0.15.0.

## Notes on versioning and semantic versioning guarantees

Unlike sourmash,
[which provides guarantees that command-line options and outputs will not change within minor versions](https://sourmash.readthedocs.io/en/latest/support.html#versioning-and-stability-of-features-and-apis),
we make no guarantees of stability within the branchwater plugin. This
is because the branchwater plugin is intended to move fast and
occasionally break things.

Eventually we expect to provide all of the branchwater plugin's functionality within the sourmash package, at which time the sourmash guarantees will apply!

However, we do not expect command line options and output file formats
to change quickly.

We will also endeavor to avoid changing column names in CSV output, although, we may change the _order_ of column names on occasion. Please use the column headers (column names) to select specific columns.

## Notes on performance, concurrency and efficiency

The branchwater plugin is focused on power use cases and supporting
tradeoffs between disk I/O, memory, and CPU. That means things get
complicated :).

Each branchwater command does things somewhat differently, with
implications for CPU and disk load; moreover, the performance will
vary depending on the database format. You can measure threading
efficiency with `/usr/bin/time -v` on Linux systems, and disk load by
number of complaints received when running.  Your best bet is to
[just ask the team for suggestions](https://github.com/dib-lab/sourmash/issues),
but you should lightly skim the docs below and play with some small data
sets, too!

However, perhaps the single biggest boost to efficiency can come from
choosing different (higher) scaled parameters: you can set the scaled
parameter higher if you're only looking for large overlaps. For example,
if you're looking for overlaps > 50kb, you can set scaled to 10000
(`-s 10000`) and most commands should get at least 10x faster.

---

`manysketch` loads one sequence file from disk per thread and sketches
it using all signature params simultaneously.

`manysearch` loads all the queries at the beginning, and then loads
one database sketch from disk per thread. The
compute-per-database-sketch is dominated by I/O. So your number of
threads should be chosen with care for disk load. We typically limit
it to `-c 32` for shared disks. <!-- We suggest using a manifest CSV
file for the database sketches. CTB -->

`multisearch` loads all the queries and database sketches once, at the
beginning, and then uses multithreading to search across all matching
sequences. For large databases it is extremely efficient at using all
available cores. So 128 threads or more should work fine!
Zipfiles <!-- and manifests CTB --> should work well.

`pairwise` acts just like `multisearch`, but only loads one file (and
then does all comparisons between all pairs within that file).

Like `multisearch` and `pairwise`, `fastgather` loads everything at
the beginning, and then uses multithreading to search across all
matching sequences. For large databases it is extremely efficient at
using all available cores. So 128 threads or more should work fine! We
suggest using zipfiles <!-- or manifests CTB --> for the database.

`fastmultigather` loads the entire database once, and then loads one
query from disk per thread. The compute-per-query can be significant,
though, so multithreading efficiency here is less dependent on I/O and
the disk is less likely to be saturated with many threads. We suggest
limiting threads to between 32 and 64 to decrease shared disk load.

`cluster` loads the entire file multithreaded, and then populates the
graph sequentially.

## Appendix 1 - `index` to create a low-memory index

The command `sourmash scripts index` makes an on-disk inverted index
for low memory fast search. Indexing takes a while, but then search
takes fewer resources.

Currently only `fastmultigather` and `manysearch` can use this kind of index.

`fastmultigather` with this index produces a complete set of `sourmash gather` columns.

We suggest using the extension `.rocksdb` for these databases, as we
use [RocksDB](https://rocksdb.org/) for the underlying database storage
mechanism.
