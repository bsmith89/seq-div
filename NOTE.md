---
title: Methanogen diversity of the human gut
author: Byron J. Smith
...
<!-- Top Matter {{{-->

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

<!-- /Top Matter }}}-->

# Notebook #
## Clone Library Sequence Analysis ##
(date: 2015-05-04)

(date: 2015-05-05)

(date: 2015-05-06)

-  Recieved sequences from core.
-  Analyzed diversity and found that 3 clades of methanogens were present in
   my clones.
-  Individuals, in general, had just one clade.
-  I even found some people with Methanomasillococcus, which is really cool.

![This tree shows the outstanding clade, and the segregation by individual.](static/2015-05-07_tree.png)

(date: 2015-05-07)

-  Based on the cleaned up alignments, the ATA C-terminus of the consensus
   protein alignment does not match with (any of?) the reference sequences.
   This sequence appears to overhang on the 5' ends of the luton amplicons I
   pull out from my clone sequences.

![See the clones (at the top) have an aligned 5' overhang that doesn't match _any_ of the references](static/2015-05-07_alignment.png)


# Appendices #
## Data Sources ##
### [`raw/*.ab1`] ###
Raw traces were downloaded from the UM sequencing core webserver.

These are binary files in the ABI format.

Eventually, this file should be publicly available online.  For now, it is
on the sequencing core servers.

```bash
# This file has a recipe in Makefile.
make download-traces
```

## Data Semantics ##
### [`res/*.primersearch-luton.tsv`] ###
[Explanation for column meanings.]

[Any notes regarding individual observations.]
