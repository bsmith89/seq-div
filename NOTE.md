---
title: Methanogen diversity of the human gut
author: Byron J. Smith
...

# Introduction #
Sequence based analysis of the composition of the mcrA gene population.

[A larger passage detailing the motivation and direction of the project.]

## Requirements and Environment ##
### Applications ###
-  [Python](http://www.python.org/) (3.4+)
-  GNU Make
-  BASH
-  phred
-  Muscle
-  FastTree

# Notebook #
## Clone Library Sequence Analysis ##
(date: 2015-05-04 through 2015-05-06)

-  Recieved sequences from core.
-  Analyzed diversity and found that 3 clades of methanogens were present in
   my clones.
-  Individuals, in general, had just one clade.
-  I even found some people with Methanomasillococcus, which is really cool.

(date:[Continuation Date])

# Appendices #
## Data Sources ##
### [`raw/*.ab1`] ###
Raw traces were downloaded from the UM sequencing core webserver.

These are binary files in the ABI format.

Eventually, this file should be publicly available online.  For now, it is
on the sequencing core servers.

```bash
# This file has a recipe in Makefile.
make raw/example-raw-file.txt
```

## Data Semantics ##
### [`res/*.primersearch-luton.tsv`] ###
[Explanation for column meanings.]

[Any notes regarding individual observations.]
