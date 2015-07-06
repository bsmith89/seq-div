# Preface {{{1
# User help message {{{2
define HELP_MSG

================================
 Analysis Makefile Documentation
================================

SYNOPSIS
    Run project operations using make commands.

TARGETS
    all
        Generate all. (By default this includes all figures,
        results, and documentation based on user-defined recipes.)

    docs
        Compile markdown files (e.g. NOTE.md, TEMPLATE.md) into HTML (uses
        Pandoc).

    figs
        Carry out the pipeline to ultimately generate figures of the results.

    res
        Carry out the pipeline to ultimately generate quantitative results
        files.

    help
        Show this help message.

    init
        Initialize the project:
            (1) submodules
            (2) venv
            (3) python-reqs
            (3) data-dirs
            (4) configure git to automatically clean IPython notebooks;
            (5) OPTIONAL: remove the 'origin' git remote
            (6) OPTIONAL: squash the commit history into a single
                'Initial commit';
            (7) create `.git/.initialized` to indicate that these steps are
                completed.

    submodules
        Initialize and update all git submodules (see `.gitmodules`).

    venv
        Create the virtualenv if absent.

    python-reqs
        Install all python requirements from requirements.txt and
        all <SUBMODULE>/requirements.txt to the venv.

    data-dirs
        Create all data directories listed in $${DATA_DIRS}
        Default: raw/ seq/ tre/ img/ fig/ res/

EXAMPLES
    make init  # Initialize the project.
    make all   # Carry out all defined steps in the project.

GNU MAKE HELP:


endef
export HELP_MSG


# ========================
#  Standard Configuration {{{2
# ========================
# One failing step in a recipe causes the whole recipe to fail.
.POSIX:

# Don't delete intermediate files.
.SECONDARY:

# Delete targets if there is an error while executing a rule.
.DELETE_ON_ERROR:

# The target `all` needs to be the first one defined (besides special
# targets) in order for it to be made on running `make` without a target.
.PHONY: all
all:

HELP_TRGTS = help h HELP Help
.PHONY: ${HELP_TRGTS}
${HELP_TRGTS}:
	@echo "$$HELP_MSG" "$$(${MAKE} -h)" | less

# All recipes are run as though they are within the virtualenv.
# WARNING: This may cause difficult to debug problems.
VENV = ./venv
export VIRTUAL_ENV = $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}

# TODO: Include a tmp/ dir?  Use it for what?
DATA_DIRS += etc/ ipynb/ raw/ meta/ res/ fig/

# Use this file to include sensitive data that shouldn't be version controlled.
# Others forking this project will need to create their own local.mk.
# If local.mk is vital and you would like the user to be alerted to its
# absence, remove the preceeding '-'.
include local.mk

# ====================
#  User Configuration {{{1
# ====================

# Major targets {{{2
.PHONY: docs figs res
docs:
figs:
res: tre/mcra-both2.luton-ampli.qtrim.gb.nucl.nwk \
     tre/mcra-both2.luton-ampli.qtrim.gb.prot.nwk \
     tre/mcra-refs.gb.prot.nwk \
     tre/rrs-refs.gb.nucl.nwk

# What files are generated on `make all`?
all: docs figs res

# Compute environment {{{2
# Name, and directory, of the python virtual environment:
VENV = ./venv
# All recipes are run as though they are within the virtualenv.
# WARNING: This may cause difficult to debug problems.
# To deactivate, thereby running all recipes from the global python
# environment, comment out the following line:
export VIRTUAL_ENV = $(abspath ${VENV})

# Use the following line to add to the PATH of all recipes.
# WARNING: These executibles will not necessarily be available in the same
# way from the command line, so you may get difficult to debug problems.
export PATH := ${VIRTUAL_ENV}/bin:${PATH}
# TODO: Deal with virtualenvs in a more transparent way.

# Cleanup settings {{{2
# Use the following line to add files and directories to be deleted on `make clean`:
CLEANUP += res/* seq/* tre/* meta/*  # Remove intermediate results
CLEANUP += raw/unarchive.mk raw/all.fastq.mk  # Remove sub-makefiles

# Initialization settings {{{2
# What directories to generate on `make data-dirs`.
# By default, already includes etc/ ipynb/ raw/ meta/ res/ fig/
DATA_DIRS += seq/ tre/ raw/ab1

# ==============
#  Data {{{1
# ==============
# User defined recipes for cleaning up and initially parsing data.
# e.g. Slicing out columns, combining data sources, alignment, generating
# phylogenies, etc.

# Archives {{{2
# Dropbox downloads {{{3

define GET_FROM_MY_DROPBOX
wget --no-clobber --no-verbose --directory-prefix=${@D} \
	https://dl.dropboxusercontent.com/u/$${DROPBOX_UID}/Data/${@F}
endef
# TODO: Can I count on all raw/*.tgz files being in my dropbox?
raw/%.tgz:
	${GET_FROM_MY_DROPBOX}

raw/mcra.published.fn:
	${GET_FROM_MY_DROPBOX}

# Unpacking and copying {{{3
# Automatically generate mapping between .ab1 sequences and their archive.
raw/ab1.mk: etc/mcra-clones.meta.tsv
	awk 'NR > 1 {print "raw/ab1/" $$1 ".ab1: raw/" $$9}' < $^ > $@
include raw/ab1.mk

# And now the actual recipe
# TODO: Since we only take the first pre-requisite as the file source,
# problems in which multiple preqs are set may be difficult to debug.
raw/ab1/%.ab1:
	cp $</${@F} $@

# Automatically generate rules for unpacking archived data files.
raw/unarchive.mk: etc/mcra-clones.meta.tsv
	awk 'NR > 1 && !seen[$$9]++ {print "raw/" $$9 ": raw/" $$9 ".tgz\n\ttar -C raw/ -xzf $$^\n\ttouch $$@"}' < $^ > $@
include ./raw/unarchive.mk

# Contruct FASTQs {{{2
# Process ab1's {{{3
raw/ab1/%.ab1.seq raw/ab1/%.ab1.qual: raw/ab1/%.ab1
	phred $< -qd raw/ab1 -sd raw/ab1 -raw $* >/dev/null

raw/ab1/%.fastq: bin/make_fastq.py raw/ab1/%.ab1.seq raw/ab1/%.ab1.qual
	$(word 1,$^) $(word 2,$^) $(word 3,$^) > $@

# Concatenate clone FASTQ's {{{3
# All pre-requisites for raw/mcra-clones.all.fastq
raw/all.fastq.mk: etc/mcra-clones.meta.tsv
	awk 'NR>1 {print "raw/mcra-clones.all.fastq: raw/ab1/" $$1 ".fastq"}' $^ > $@
include raw/all.fastq.mk

# And now the actual recipe
raw/mcra-clones.all.fastq:
	cat $^ > $@

# Process clone metadata {{{2
meta/mcra-clones.annot.tsv: etc/mcra-clones.meta.tsv
	awk '{print $$2 "." $$3 "." $$4 "." $$5 "." $$6 "." $$7 "\t" $$0 }' $^ > $@

meta/mcra-clones.names.tsv: meta/mcra-clones.annot.tsv
	awk 'NR > 1 {print $$2 "\t" $$1}' $^ > $@

meta/mcra-clones.list: meta/mcra-clones.annot.tsv
	awk 'NR > 1 && $$9 != "True" {print $$1}' $^ > $@
# Include only the clones which were amplified with a particular set of primers
# and which are not suspect.

# Produce a curated set of clones {{{2
seq/mcra-clones.fastq: raw/mcra-clones.all.fastq \
                             bin/utils/rename_seqs.py meta/mcra-clones.names.tsv \
                             bin/utils/fetch_seqs.py meta/mcra-clones.list
	cat $(word 1,$^) \
        | $(word 2,$^) -f fastq -t fastq $(word 3,$^) \
        | $(word 4,$^) -f fastq -t fastq $(word 5,$^) > $@

# Reference sequences {{{2
# Rename reference sequences and remove those which have been a priori deemed
# suspicious.

# Process ref metadata {{{3
# TODO: Why can't I import these annotations into FigTree?
meta/%-refs.annot.tsv: etc/refs.meta.tsv
	awk '$$3 == "$*" || NR == 1 {print $$2 "." $$1 "\t" $$0}' $^ > $@

meta/%-refs.names.tsv: meta/%-refs.annot.tsv
	awk 'NR > 1 {print $$2 "\t" $$1}' $^ > $@

meta/%-refs.suspect.list: meta/%-refs.annot.tsv
	awk '$$4 == "$*" && $$6 == "True" {print $$1}' $^ > $@

# Make curated mcrA reference sets {{{3
# De-redundant-ed references (mcra-refs.names.tsv isn't exhaustive).
seq/mcra-refs.fn: bin/utils/rename_seqs.py meta/mcra-refs.names.tsv \
                  raw/mcra.published.fn \
                  bin/utils/drop_seqs.py meta/mcra-refs.suspect.list
	$(word 1,$^) $(word 2,$^) $(word 3,$^) \
        | $(word 4,$^) $(word 5,$^) > $@

# Very reduced subset of reference sequences which guide interpretation of
# my clone libraries.
meta/%-refs.reduced.list: meta/%-refs.annot.tsv
	awk '$$4 == "$*" && $$5 == "True" {print $$1}' $^ > $@

seq/mcra-refs2.fn: bin/utils/rename_seqs.py meta/mcra-refs.names.tsv \
                   raw/mcra.published.fn \
                   bin/utils/fetch_seqs.py meta/mcra-refs.reduced.list
	$(word 1,$^) $(word 2,$^) $(word 3,$^) \
        | $(word 4,$^) $(word 5,$^) > $@

# Curated 16S references {{{3
# rrs-refs3 is a synonym for rrs-refs
seq/rrs-refs.fn: bin/utils/rename_seqs.py meta/rrs-refs.names.tsv \
                 raw/rrs.published.fn \
                 bin/utils/drop_seqs.py meta/rrs-refs.suspect.list
	$(word 1,$^) $(word 2,$^) $(word 3,$^) \
        | $(word 4,$^) $(word 5,$^) > $@


# 16S Reference alignment/hmm {{{3
raw/Silva.seed_v119.tgz:
	curl -o $@ http://www.mothur.org/w/images/5/56/Silva.seed_v119.tgz

raw/silva.seed_v119.align: raw/Silva.seed_v119.tgz
	tar -C ${@D} -xzf $^ ${@F}
	touch $@

res/rrs.hmm: raw/silva.seed_v119.align
	hmmbuild --dna --informat afa $@ $^


# Combine clones with references {{{3
# Combine clones which have been quality trimmed with the reference sequences.
# to make the "both" file series.
# Quality trimming of the references is not required.
seq/mcra-both.%-ampli.qtrim.fn: seq/mcra-clones.%-ampli.qtrim.fn seq/mcra-refs.%-ampli.fn
	cat $^ > $@

seq/mcra-both2.%-ampli.qtrim.fn: seq/mcra-clones.%-ampli.qtrim.fn seq/mcra-refs2.%-ampli.fn
	cat $^ > $@

seq/mcra-both2.f3r4-ampli.qfilt.afn: seq/mcra-clones.f3r4-ampli.qfilt.afn seq/mcra-refs2.f3r4-ampli.afn
	cat $^ > $@

# Remove uniformative sequence {{{2
# Excise amplicon {{{3
# Search for primer hits:
res/%.psearch.out: seq/%.fastq etc/primers.tsv
	primersearch -seqall $(word 1,$^) -sformat fastq -infile $(word 2,$^) -mismatchpercent 30 -outfile $@
# I use 30% mismatchpercent so that I can include many hits before
# narrowing by other criteria.

res/%.psearch.out: seq/%.fn etc/primers.tsv
	primersearch -seqall $(word 1,$^) -infile $(word 2,$^) -mismatchpercent 20 -outfile $@
# I use 20% mismatchpercent since quality isn't a concern coming from
# as FASTA (not FASTQ) file.

res/%.psearch.tsv: bin/parse_psearch.py etc/primers.tsv res/%.psearch.out
	$(word 1,$^) $(word 2,$^) $(word 3,$^) > $@

# Re-orient the sequences to match the forward-reverse in etc/*.primers
# and trim to within the primers
# TODO: Should this be systematized somehow?  DRY.
seq/%.luton-ampli.fastq: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fastq
	$(word 1,$^) --primer-set luton --max-mismatch 5 --trim-primers --drop -f fastq -t fastq $(word 2,$^) $(word 3,$^) > $@

seq/%.luton-ampli.fn: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fn
	$(word 1,$^) --primer-set luton --max-mismatch 5 --trim-primers --drop $(word 2,$^) $(word 3,$^) > $@

seq/%.f3r4-ampli.fastq: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fastq
	$(word 1,$^) --primer-set f3r4 --max-mismatch 1 --trim-primers --drop -f fastq -t fastq $(word 2,$^) $(word 3,$^) > $@

# Based on the reference alignment, there are no gaps in this region,
# so the amplicon should be already aligned, assuming it's the right length
seq/%.f3r4-ampli.afn: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fn bin/drop_missized.py
	$(word 1,$^) --primer-set f3r4 --max-mismatch 1 --trim-primers --drop $(word 2,$^) $(word 3,$^) \
        | $(word 4,$^) \
        > $@

seq/%.f3r4-ampli.fastq: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fastq
	$(word 1,$^) --primer-set f3r4 --max-mismatch 0 --trim-primers --drop -f fastq -t fastq $(word 2,$^) $(word 3,$^) > $@

# Based on the reference alignment, there are no gaps in this region,
# so the amplicon should be already aligned
seq/%.f3r4-ampli.qfilt.afn: ./bin/qtrim_reads.py seq/%.f3r4-ampli.fastq bin/drop_missized.py
	$(word 1,$^) --just-filter $(word 2,$^) \
        | $(word 3,$^) \
        > $@

# Quality trim {{{3
seq/%.qtrim.fn: bin/qtrim_reads.py seq/%.fastq
	$^ --drop > $@

# Convert filetypes {{{2
seq/%.fn: bin/utils/convert.py seq/%.fastq
	$^ -f fastq > $@

# Translate {{{2
seq/%.frame.fn: ./bin/infer_frame.py seq/%.fn
	$^ > $@

seq/%.fa: bin/utils/translate.py seq/%.frame.fn
	$^ > $@

#  Analysis {{{1
# =======================
# User defined recipes for analyzing the data.
# e.g. Calculating means, distributions, correlations, fitting models, etc.
# Basically anything that *could* go into the paper as a table.

# Align {{{2
# refinement can't handle hmmalign output where the sequence just doesn't
# match the model at all.
seq/mcra-%.afa: etc/mcra.fungene.hmm seq/mcra-%.fa bin/utils/convert.py
	hmmalign --amino --informat FASTA $(word 1,$^) $(word 2,$^) \
        | $(word 3,$^) --in-fmt stockholm --out-fmt fasta \
        | muscle -quiet -refine \
        > $@

# TODO: Compare this to SINA
seq/rrs-%.afn: res/rrs.hmm seq/rrs-%.fn bin/utils/convert.py
	hmmalign --dna --informat FASTA $(word 1,$^) $(word 2,$^) \
        | $(word 3,$^) --in-fmt stockholm \
        | muscle -quiet -refine \
        > $@

seq/%.afa: seq/%.fa
	muscle -quiet < $^ > $@

seq/%.afn: bin/utils/codonalign.py seq/%.afa seq/%.frame.fn
	$^ > $@

# Gblocks {{{2
# gblocks can't handle sequence names longer than 72 characters
# gblocks always throws an error code of 1
seq/%.gb.afn res/%.gb.html: seq/%.afn
	Gblocks $^ -t=c -p=y -v=150 || [ $$? == 1 ]
	mv $^-gb seq/$*.gb.afn
	mv $^-gb.htm res/$*.nucl.gb.html

seq/%.gb.afa: seq/%.afa
	Gblocks $^ -t=p -p=y -v=150 || [ $$? == 1 ]
	mv $^-gb seq/$*.gb.afa
	mv $^-gb.htm res/$*.prot.gb.html

# =======================
# Trees {{{2
tre/%.nucl.nwk: seq/%.afn
	FastTree -quiet -nt $^ > $@

tre/%.prot.nwk: seq/%.afa
	FastTree -quiet $^ > $@

# ==================
#  Graphing {{{1
# ==================
# User defined recipes for plotting figures.  These should use
# the targets of analysis recipes above as their prerequisites.


# =======================
#  Documentation {{{1
# =======================
ALL_DOCS = TEMPLATE NOTE
ALL_DOCS_HTML = $(addsuffix .html,${ALL_DOCS})
MATHJAX = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

%.html: %.md
	pandoc -f markdown -t html5 -s --highlight-style pygments \
        --mathjax=${MATHJAX} --toc --toc-depth=4 --css static/main.css $^ > $@

docs: ${ALL_DOCS_HTML} fig/Makefile.reduced.png

# Visualize makefile with cytoscape.
# requires:
# https://bitbucket.org/jpbarrette/makegrapher
# details at: https://code.google.com/p/makegrapher
res/Makefile.complete: Makefile
	${MAKE} --makefile=$^ -npr > $@

res/Makefile.dot: res/Makefile.complete
	make_grapher.py -T $^ -o $@ >/dev/null

res/Makefile.reduced.dot: bin/clean_makefile_graph.py res/Makefile.dot
	$(word 1,$^) -d '^raw/ab1' -d '^bin/utils/' -d '^\.' -d '\.git' \
                 -d '(submodules|venv|python-reqs|init)' \
                 -k '^raw/mcra' -k '^(all|res|figs|docs|Makefile)$$' \
                 $(word 2,$^) > $@

fig/Makefile.reduced.png: res/Makefile.reduced.dot
	dot -Tpng -Grankdir=BT -Nshape=plaintext < $^ > $@

# =================
#  Cleanup {{{1
# =================
.PHONY: clean
clean:
	rm -rf ${CLEANUP}

# ========================
#  Initialization Recipes {{{1
# ========================
.PHONY: init
init: .git/.initialized
.git/.initialized:
	@${MAKE} submodules
	@[ "${VENV}" ] && ${MAKE} ${VENV}
	@${MAKE} python-reqs
	@${MAKE} data-dirs
	@${MAKE} .link-readme
	@${MAKE} .ipynb-filter-config
	-@${MAKE} .git-mangle
	touch $@

define VENV_ACTIVATE_MSG

A python3 virtual environment has been made in `${VENV}`.

Python called from recipes in `Makefile` will automatically use this virtual
environment.  To activate ${VENV} for the command-line, however,
run `source ${VENV}/bin/activate`.

endef
export VENV_ACTIVATE_MSG

${VENV}:
	python3 -m venv $@
	@echo "$$VENV_ACTIVATE_MSG"

# Git Submodules:
SUBMODULE_DIRS := $(shell git submodule | sed 's:^ ::' | cut -d" " -f2)
SUBMODULES = $(addsuffix /.git,${SUBMODULE_DIRS})

.PHONY: submodule python-reqs data-dirs
submodules: ${SUBMODULES}
${SUBMODULES}: .gitmodules
	git submodule update --init --recursive ${@D}

SUBMODULE_PIP_REQS = $(wildcard $(addsuffix /requirements.txt,${SUBMODULE_DIRS}))
PIP_REQS = requirements.txt ${SUBMODULE_PIP_REQS}

python-reqs: | ${VENV}
	for req_file in ${PIP_REQS} ; do \
        pip install --upgrade --no-deps -r $$req_file ; \
        pip install -r $$req_file ; \
    done

data-dirs:
	mkdir -p ${DATA_DIRS}

.PHONY: .link-readme .confirm-git-mangle \
        .git-mangle .ipynb-filter-config
.link-readme:
	unlink README.md
	ln -s NOTE.md README.md

define INIT_OPTS_MSG

You are about to remove the remote repository labeled 'origin' and squash the
commit history into a single commit.  This procedure makes sense for
initializing a new project from a template project, where the history of the
template is unimportant and you do not want to push or pull changes to the
remote repository from which you cloned the template.

Alternatively, if you are initializing a previously started project, you most
likely do not want to lose the commit history, and you may want to push or pull
changes from the remote.  In that case, respond with something other than "y"
or "Y" to the following prompt.  `make` will thow an error, but initialization
will work as intended.

endef
export INIT_OPTS_MSG

.confirm-git-mangle:
	@echo "$$INIT_OPTS_MSG"
	@read -rp "Are you sure you want to remove the remote and squash the commit history? [y/N]: " MANGLE ; \
	[ $$MANGLE == "y" ] || [ $$MANGLE == "Y" ]

.git-mangle: .confirm-git-mangle
	git remote remove origin
	git branch -m master
	git reset --soft $$(git rev-list --max-parents=0 HEAD)
	git add -A
	git commit --amend -m "Initial commit."

%/requirements.txt: %/.git

# IPython Notebook Output Filter Configuration
# TODO: Should I require `python-reqs`?
bin/utils/ipynb_output_filter.py: bin/utils/.git

.ipynb-filter-config: bin/utils/ipynb_output_filter.py
	git config --local filter.dropoutput_ipynb.clean bin/utils/ipynb_output_filter.py
	git config --local filter.dropoutput_ipynb.smudge cat
