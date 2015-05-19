---
title: Methanogen diversity of the human gut
author: Byron J. Smith
...
<!-- Top Matter {{{-->

# Introduction #
Analysis of sequence diversity of enzymes in various hydrogenotrophic
pathways in the human gut.

While it is thought that competition for hydrogen plays a large role in
determining the prevalence of methanogens in the human gut, it is difficult
to measure the selective pressures acting on hydrogenotrophic taxa.
Here I propose to use sequence data to detect signals of selection.

Sequence based analysis of the composition of the mcrA gene population.

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
## Public Sequence Analysis ##
Publicly available mcrA sequences were analyzed.

### _mcrA_ Sequences Retrieval ###
(date:2015-04-06)

Genomes were selected for two sets of methanogens on IMG.
First, a "Find Functions" search for 2.8.4.1 was carried out, yielding 115
genomes with _mcrA_.
These were added to the genome cart.
Then, a search on these genomes for "2.8.4.1" yielded 431 genes, including the
alpha, beta, and gamma subunits of methyl-coenzyme M reductase.
A subsequent narrowing by entries which had "alpha subunit" in "Gene Product
Name" field, yielded 143 sequences.
As indicated by the number of genes relative to the number of genomes, a number
of genomes have two copies of the gene.
These 143 sequences were downloaded in FASTA format as
`raw/mcra.methanogens.fn`
(Actually, they were copied from the output page with spaces in between FASTA
formatted entries and these spaces were removed with `sed '/^$/d'`.)

From these 115 genomes, a search on pfam domain pfam08423 yielded 227 genes,
for both RadA and RadB.
I believe that _radA_ is the archael homologue to _recA_.
Narrowing by entries to those with two domains,
"HHH\_5(pfam14520), Rad51(pfam08423)", in the "Domains" field, yielded 113
genes.
These 113 _radA_ sequences with both domains were downloaded as a FASTA file,
`raw/rada.methanogens.fn`.

From these 113 _radA_ sequences, the 20 Gene IDs from genomes named
"Methanobrevibacter smithii TS\*" (where "\*" was replaced with a two or three
digit number followed by a letter) were collected.
These genomes correspond with the 20 genomes described in Hansen et al. 2011
of isolates collected from human feces.
Similarly, the 21 _mcrA_ sequences from these 20 genomes were also collected.

### _mcrA_ Analysis ###
(date:2015-04-07)

Nucleotide a protein alignments were made for these sequences by translating to
amino-acids, aligning these with `muscle`, and then "backaligning" the
nucleotide sequence to match the protein alignment.
While this may have _some_ shortcomings relative to a true codon alignment,
I can't imagine they are more than minimal.
Trees were then constructed for both nucleotide and protein alignments using
`fasttree`.
The resulting branch lengths could then be output, although it is worth
mentioning that branch lengths are suspect with this FastTree.
([See here.](http://www.microbesonline.org/fasttree/#BranchLen))

The _radA_ trees have one long, suspicious branch, and the two protein
sequences at its tip not blast particularly close to any other _radA_.

![Suspicious sequences on protein tree.](static/2015-04-07.rada_prot_tree.png)

I have removed these sequences from analysis, but the _radA_ protein to nucleotide
branch length ratio is still about the same as the _mcrA_ ratio for all
methanogens.
In fact, it's suprising how close to the same total branch lengths for the two
trees are.

I'm interested in checking out other marker genes.
[This](https://phylosift.wordpress.com/tutorials/scripts-markers/old-pmprok-marker-gene-names/)
seems like a potentially useful list.
Marker genes aren't necessarily the best place to look, though.
I'm curious about the range of branch length ratios (BLRs) that I might find
in the methanogen or Hansen sets.

## Primer Selection for Methanogen Phylogenetic Marker Libraries ##

(date: 2015-04-08)

While _mcrA_ would be an ideal marker gene to explore methanogen diversity
in the human gut, it's possible that primer design will be difficult for this
gene.
Alternatively, I could use archael 16S sequences but these potentially
have too little variation to observe strain-level diversity.
Genomes from Hansen et al. have very little diversity in _mcrA_, _mcrB_,
_mcrG_, etc. but this could be due to the process of isolation.
If I want to say something about strain diversity in (non-)methane producing
human subjects, I will want to use culture independent approaches.
Jessica Seiber has a set of primers which will amplify a region of 450-500bp.
This is probably too big for Illumina sequencing, but would work for a clone
library.

These primers are from Luton 2002:

| Forward: `GGTGGTGTMGGATTCACACARTAYGCWACAGC`
| Reverse: `TTCATTGCRTAGTTWGGRTAGTT`

I took a stab at designing my own primers. See `ipynb/primers.ipynb`.
The first result was this figure:

![Inferences about potential primers from an _mcrA_ alignment.](static/2015-04-08.primer_design.png)

Here we see that there are regions at approximately 550, 1100, 1500, and 1550
where a 15 base primer could have a degeneracy of less than 100 (ignoring
sequence variants at a position with only one representative).
An attempt to design primers at these positions was found to be very difficult.
As a potentially easier task, primers were designed for the two
branches separated by the the deepest split in
the phylogenetic tree with 82 and 61 taxa represented, respectively.

![Finding primer sites in two clades of _mcrA_.](static/2015-04-08.primer_design2.png)

This was done by hand with an attempt to keep degeneracy below 100, primer
lengths of ~20 bases, and ignoring sequence variants at a position with fewer
than 5 representatives.
The result were the following sequences at primer locations
(on the coding strand, so the reverse primer has to be reverse complemented):

Split A (82 taxa)
| 561 Forward: `CARGARCAYATGGTDGAR`
| 1592 Reverse: `CCWAACTAYCCWAACTAYGCHATGAAYGT`

Split B (61 taxa)
| 372 Forward: `GAYGAYCTBCACTWYGTSAACAA`
| 558 Forward: `GTBCAGGARMWSATGG`  <- Probably too short
| 1583 Reverse: `GGHSCVAACTAYCCSAACTAYGC`

I also spent some time automating the primer design process myself, and
had Tom S. recommend the software, "CODEHOP".

(date: 2015-04-10)

Potential primers can be analyzed with
[OligoAnalyzer 3.1 from IDT](https://www.idtdna.com/calc/analyzer).
Other resources can be found
[here](http://bitesizebio.com/18992/a-primer-for-designing-degenerate-primers/),
and [here](http://www.lifetechnologies.com/us/en/home/products-and-services/product-types/primers-oligos-nucleotides/invitrogen-custom-dna-oligos/primer-design-tools.html).

Methanogenesis is defined for KEGG Module M00567.

## Clone Library Sequence Analysis ##
(date: 2015-04-28)

A single test clone came back with sequence.
When I blast the sequence, the best hit stretches from base 17 to 485.
The matching sequence goes all the way from the Luton Forward to Luton Reverse
primers, and is therefore a full length sequence of the same clone that I
sequenced.
While there are some errors near the 3' end of my sequence (a single InDel in
the primer sequence), this means that I got back the longest length I could
have hoped.

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

(date: 2015-05-19)

Sequence "2481783.M13REV" which I have listed as a re-sequence of U003.2.1.C03
appears to actually be a re-sequence of "U003.2.1.C02".
Which, if I'm not mistaken, is a good thing, since it has NO sequences, at
least within the g-blocks, and yet _does_ differ from many other sequences,
e.g. U011.3.1.B01 at two positions.

I'll correct the meta-data, but I have to correct this in my notebook.
The two mutations are non-synonymous.


### Things to present at lab meeting ###

1.  Two species of mcrA possessing organisms in the human gut
2.  Two clades of _smithii_.
3.  Clades segregate into individuals
4.  Some nucleotide diversity in _M. smithii_
    A.  Seems biased towards third position nucleotides
    B.  But sample size pretty low and may be erroneous
5.  Caveats
    A.  Single PCR per sample
    B.  Low sample size
         -  Subjects
         -  Clones
    C.  Sequencing error
    D.  Chimeras
    E.  Primer bias/limitation

### Possible experimental design ###

 -  Design primers that will work for full-overlap paired ends Illumina
 -  Two independent PCR runs for each sample
     -  Only consider sequences which appear in both runs

(date: 2015-05-12)
TODO: I need to see if the separation of M. smithii into two clades is also
visible in 16S sequence.


# Appendices #
## Data Sources ##
### `raw/<repository>/*.ab1` ###
Raw traces were downloaded from the UM sequencing core webserver individually.

These are binary files in the ABI format.

This file is be publicly available online at:
https://dl.dropboxusercontent.com/u/\<UID\>/\<PATH_TO_FILE_IN_PUBLIC_FOLDER\>

These files can be re-downloaded automatically based on a rule in `Makefile`.

### `raw/mcra.refs.fn` ###
All mcrA sequences were downloaded from IMG.
Sequences were obtained from IMG/ER last updated 2015-04-02.
The page footer says 'Version 4.510 Oct 2014'
Sequences were found using Find Functions -> Enzyme (KEGG) and searching
for 2.8.4.1
(methyl-coenzyme M reductase a.k.a Coenzyme-B sulfoethylthiotransferase).
115 genomes were found which had this function.
Searching these genomes for 2.8.4.1 yielded 431 genes, some of which were
subunits beta and gamma, rather than just alpha.
143 sequences for subunit beta were downloaded.

Files are FASTA formatted with labels starting with a unique (to IMG) gene
ID and followed by a full description of the sequenced organism and the name
of the gene.

This file can be re-generated from IMG by following the instructions above.

## Data Semantics ##
### [example] ###
[Explanation for column meanings.]

[Any notes regarding individual observations.]
