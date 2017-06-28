# circtools

## What is it for?

`circtools` R package is intended to simplify the visualisation of circular 
and linear transcripts and to design primers which span through the 
splice junctions in order to catch circs in your experiments.

## What do I need to provide?

#### Circs

First of all, the genomic coordinates of splice junctions for circular RNA
transcripts.

#### Annotation 

The implementation leans on Bioconductor packages for working with 
*ENSEMBL* annotation. It is easier if you are going to use `circtools` for
organisms for which `ensembldb` packages are available,
for example, `EnsDb.Hsapiens.v86` or `EnsDb.Mmusculus.v79`.
Otherwise, it is always possible to create your own `ensembldb` object.

[How to create my own annotation](https://bioconductor.org/packages/devel/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#1022_from_a_gtf_or_gff_file) is well documented on
the `ensembldb` Bioconductor page.

## How to use

[Documentation page](../../doc/README.md)

## FAQ

### I don't have permissions to install a package on a server/machine.

You can install the package to your home directory.
You need to create/edit a file `.Renviron` at your home dir to tell R where 
to keep the installed packages.

In the file you need to add (in case you want to keep 
all at `YOUR_HOME/R/library`):

```
## Linux
R_LIBS=~/R/library

## Windows
R_LIBS=C:/R/library
```
The home folder can be found from R itself irregardless of the platform: 

`path.expand("~")`

