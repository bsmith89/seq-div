# Preface {{{1
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
#  Standard Configuration
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
VENV = ./venv
export VIRTUAL_ENV = $(abspath ${VENV})
export PATH := ${VIRTUAL_ENV}/bin:${PATH}

DATA_DIRS = etc/ ipynb/ raw/ meta/ res/ fig/

# Use this file to include sensitive data that shouldn't be version controlled.
-include local.mk

# ====================}}}
#  User Configuration {{{0
# ====================
# Use the following line to add project directories with executibles
# to the `make` recipe path:
# export PATH := <BIN-DIR-A>:<BIN-DIR-B>:${PATH}

# Use the following line to add files to be deleted on `make clean`:
CLEANUP = ${ALL_DOCS_HTML}

# What directories to generate on `make data-dirs`.
# By default, already includes etc/ ipynb/ raw/ meta/ res/ fig/
DATA_DIRS += seq/ tre/

# Add sub-targets (prerequisites) to the major phony targets.
.PHONY: docs figs res
docs:
figs:
res:

# What files are generated on `make all`?
all: docs figs res

# }}}
# ==============
#  Data Recipes
# ==============
# User defined recipes for cleaning up and initially parsing data.
# e.g. Slicing out columns, combining data sources, alignment, generating
# phylogenies, etc.

# Download and extract usable data {{{1
RAW_CLONES_SEQ_NAMES := $(shell cut -f1 etc/clones.names.tsv)
RAW_CLONES_AB1 = $(patsubst %,raw/%.ab1,${RAW_CLONES_SEQ_NAMES})
REPOSITORY_URL_BASE = http://seqcore.brcf.med.umich.edu/users/schmidt/smith
${RAW_CLONES_AB1}:
	wget --user=$$SEQCORE_USER --password=$$SEQCORE_PSWD -O $@ ${REPOSITORY_URL_BASE}/${@F}

# TODO: Should I recommend using %-pattern rules whenever possible?
# or is it better to start with hard-coded rules and then switch to them
# later?
raw/%.fastq: raw/%.ab1 bin/make_fastq.py
	phred raw/$*.ab1 -qd raw/ -sd raw/ -raw $*
	bin/make_fastq.py $<.seq $<.qual > $@
	rm $<.seq $<.qual

# Concatenate all of the experimental sequences into a raw fastq file.
RAW_CLONES_FASTQ = $(patsubst %,raw/%.fastq,${RAW_CLONES_SEQ_NAMES})
seq/clones.fastq: ${RAW_CLONES_FASTQ}
	cat $^ > $@
# }}}

# Rename and remove known bad sequences {{{1
# Pre-pre-processing of the sequence files to remove subjectively identified
# suspect sequences and renaming sequences to a more meaningful scheme.
seq/%.rename.fastq: bin/utils/rename_seqs.py etc/%.names.tsv seq/%.fastq
	$(word 1,$^) --in-fmt fastq --out-fmt fastq $(word 2,$^) < $(word 3,$^) > $@

seq/%.no-susp.fastq: bin/utils/drop_seqs.py etc/%.suspect.list seq/%.fastq
	$(word 1,$^) --in-fmt fastq --out-fmt fastq $(word 2,$^) < $(word 3,$^) > $@

seq/%.rename.no-susp.fastq: bin/utils/drop_seqs.py etc/%.suspect.list seq/%.rename.fastq
	$(word 1,$^) --in-fmt fastq --out-fmt fastq $(word 2,$^) < $(word 3,$^) > $@

seq/%.fn: bin/fq2fn.py seq/%.fastq
	$^ > $@

# Rename reference sequences
seq/refs.fn: bin/utils/rename_seqs.py etc/refs.names.tsv raw/mcra.refs.fn
	$^ > $@

seq/both.ampli.qtrim.fn: seq/clones.ampli.qtrim.fn seq/refs.ampli.fn
	cat $^ > $@

# }}}

# Search for primer hits:
res/%.psearch.out: seq/%.fn etc/primers.tsv
	primersearch -seqall $(word 1,$^) -infile $(word 2,$^) -mismatchpercent 40 -outfile $@
# I use 40% mismatchpercent so that I can include many hits before
# narrowing by other criteria.

res/%.psearch.out: seq/%.fastq etc/primers.tsv
	primersearch -seqall $(word 1,$^) -sformat fastq -infile $(word 2,$^) -mismatchpercent 40 -outfile $@
# I use 40% mismatchpercent so that I can include many hits before
# narrowing by other criteria.

res/%.psearch.tsv: res/%.psearch.out bin/parse_psearch.py
	$(word 2,$^) < $< > $@

# Re-orient the sequences to match the forward-reverse in etc/*.primers
# and trim to within the primers
seq/%.ampli.fastq: ./bin/find_amplicon.py etc/primers.tsv res/clones.psearch.tsv seq/%.fastq
	$(word 1,$^) --in-fmt fastq --out-fmt fastq $(word 2,$^) $(word 3,$^) $(word 4,$^) > $@

seq/%.ampli.fn: ./bin/find_amplicon.py etc/primers.tsv res/%.psearch.tsv seq/%.fn
	$^ > $@

seq/%.qtrim.fn: ./bin/qtrim_reads.py seq/%.fastq
	$^ > $@

seq/%.fn: ./bin/fq2fn.py seq/%.fastq
	$^ > $@

seq/%.frame.fn: ./bin/infer_frame.py seq/%.fn
	$^ > $@

# Processing of clean sequences {{{2
# Translation:
seq/%.fa: bin/utils/translate.py seq/%.frame.fn
	$^ > $@

# Alignment
seq/%.afa: seq/%.fa
	cat $^ | muscle > $@

# Backalign:
# # Reordering:  *`muscle` puts sequences out of order
# seq/%.ord.afa: bin/utils/ls_ids.py bin/utils/fetch_seqs.py seq/%.afa seq/%.fn
# 	$(eval $*_TMP := $(shell mktemp))
# 	bin/utils/ls_ids.py seq/$*.fn > ${$*_TMP}
# 	bin/utils/fetch_seqs.py ${$*_TMP} seq/$*.afa > $@
# 	rm ${$*_TMP}
#
# seq/%.afn: bin/utils/backalign.py seq/%.ord.afa seq/%.fn
# 	$^ > $@

seq/%.afn: bin/utils/backalign.py seq/%.afa seq/%.frame.fn
	$^ > $@


# Trees
tre/%.nucl.nwk: seq/%.afn
	fasttree -nt $^ > $@

tre/%.prot.nwk: seq/%.afa
	fasttree < $^ > $@

# }}}

# =======================
#  Analysis Recipes
# =======================
# User defined recipes for analyzing the data.
# e.g. Calculating means, distributions, correlations, fitting models, etc.
# Basically anything that *could* go into the paper as a table.


# ==================
#  Graphing Recipes
# ==================
# User defined recipes for plotting figures.  These should use
# the targets of analysis recipes above as their prerequisites.


# =======================
#
#  Documentation Recipes {{{1
# =======================
ALL_DOCS = TEMPLATE NOTE
ALL_DOCS_HTML = $(addsuffix .html,${ALL_DOCS})
MATHJAX = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

%.html: %.md
	pandoc -f markdown -t html5 -s --highlight-style pygments --mathjax=${MATHJAX} --toc --toc-depth=4 --css static/main.css $^ > $@

docs: ${ALL_DOCS_HTML}

# =================
#  Cleanup Recipes {{{1
# =================
.PHONY: clean
clean:
	rm -f ${CLEANUP}

# ========================
#  Initialization Recipes {{{1
# ========================
.PHONY: init
init: .git/.initialized
.git/.initialized:
	@${MAKE} submodules
	@${MAKE} ${VENV}
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
