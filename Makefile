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
DATA_DIRS = etc/ ipynb/ raw/ meta/ res/ fig/

# Use this file to include sensitive data that shouldn't be version controlled.
# Others forking this project will need to create their own local.mk.
# If local.mk is vital and you would like the user to be alerted to its
# absence, remove the preceeding '-'.
include local.mk

# ====================
#  User Configuration {{{1
# ====================

# Add sub-targets (prerequisites) to the major phony targets.
.PHONY: docs figs res
docs:
figs:
res: tre/mcra-both2.luton-ampli.qtrim.gb.nucl.nwk \
     tre/mcra-both2.luton-ampli.qtrim.gb.prot.nwk

# What files are generated on `make all`?
all: docs figs res

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

# Use the following line to add files and directories to be deleted on `make clean`:
CLEANUP += res/* seq/* tre/* meta/*  # Remove intermediate results
CLEANUP += raw/unarchive.mk raw/mcra-clones.all.fastq.mk  # Remove sub-makefiles

# What directories to generate on `make data-dirs`.
# By default, already includes etc/ ipynb/ raw/ meta/ res/ fig/
DATA_DIRS += seq/ tre/

# ==============
#  Data {{{1
# ==============
# User defined recipes for cleaning up and initially parsing data.
# e.g. Slicing out columns, combining data sources, alignment, generating
# phylogenies, etc.

# Archives {{{2
# Dropbox downloads {{{3

# List of <dir>.tgz data repos in my Dropbox Public/Data/ directory.
DROPBOX_REPOS = 2015-04-29_mcrA-clones.tgz 2015-05-13_mcrA-clones.tgz 2015-06-09_mcrA-clones.tgz
define GET_FROM_DROPBOX
wget --no-clobber --directory-prefix=${@D} \
	https://dl.dropboxusercontent.com/u/$${DROPBOX_UID}/Data/${@F}
endef
$(addprefix raw/,${DROPBOX_REPOS}):
	${GET_FROM_DROPBOX}

raw/mcra.published.fn:
	${GET_FROM_DROPBOX}

# Unpacking and copying {{{3


meta/archive-contents.tsv: etc/mcra-clones.meta.tsv
	sed '1,1d' $< | awk '{print $$1 ".ab1\t" $$8}' > $@

# Automatically generate archive.mk: rules for unpacking archived data files.
raw/unarchive.mk: meta/archive-contents.tsv
	@echo Generating $@...
	@cat $^ | awk '{printf("raw/%s: raw/%s/%s\n\tcp $$^ $$@\n", $$1, $$2, $$1)}' > $@
	@for archive in $$(cut -f2 $^ | sort | uniq); do \
		printf 'raw/%s/%%.ab1 raw/%s/%%.ab1: raw/%s.tgz\n\ttar -C raw/ -xzf $$^\n\ttouch $${@D}/*\n' $$archive $$archive $$archive >> $@ ; \
	done
include ./raw/unarchive.mk

# Contruct FASTQs {{{2
# TODO: Should I recommend using %-pattern rules whenever possible?
# or is it better to start with hard-coded rules and then switch to them
# later?

raw/qual raw/seq:
	mkdir -p $@

raw/seq/%.ab1.seq raw/qual/%.ab1.qual: raw/%.ab1 | raw/qual raw/seq
	phred $< -qd ${<D}/qual -sd ${<D}/seq -raw $*

raw/%.fastq: bin/make_fastq.py raw/seq/%.ab1.seq raw/qual/%.ab1.qual
	$(word 1,$^) $(word 2,$^) $(word 3,$^) > $@

include raw/mcra-clones.all.fastq.mk
# All pre-requisites for raw/mcra-clones.all.fastq
raw/mcra-clones.all.fastq.mk: meta/mcra-clones.names.tsv
	@echo Generating $@...
	@echo > $@
	@for clone in $$(cut -f1 $^); do \
		printf 'raw/mcra-clones.all.fastq: raw/%s.fastq\n' $$clone >> $@; \
	done

raw/mcra-clones.all.fastq:
	cat $^ > $@

meta/mcra-clones.names.tsv: meta/mcra-clones.annot.tsv
	awk '{print $$2 "\t" $$1}' $^ \
		| sed '1,1d' > $@

meta/mcra-clones.annot.tsv: etc/mcra-clones.meta.tsv
	awk '{print $$2 "." $$3 "." $$4 "." $$5 "." $$6 "\t" $$0 }' $^ > $@


seq/mcra-clones.fastq: raw/mcra-clones.all.fastq \
				  bin/utils/drop_seqs.py etc/mcra-clones.suspect.list \
				  bin/utils/rename_seqs.py meta/mcra-clones.names.tsv
	cat $(word 1,$^) \
		| $(word 2,$^) -f fastq -t fastq $(word 3,$^) \
		| $(word 4,$^) -f fastq -t fastq $(word 5,$^) > $@

# Reference sequences {{{2
# Rename reference sequences and remove those which have been a priori deemed
# suspicious.

meta/mcra-refs.list: etc/mcra-refs.names.tsv
	cut -f2 $^ > $@

seq/mcra-refs.fn: bin/utils/rename_seqs.py etc/mcra-refs.names.tsv raw/mcra.published.fn \
			 bin/utils/fetch_seqs.py meta/mcra-refs.list
	$(word 1,$^) $(word 2,$^) $(word 3,$^) \
		| $(word 4,$^) $(word 5,$^) > $@

seq/mcra-refs2.fn: bin/utils/rename_seqs.py etc/mcra-refs.names.tsv raw/mcra.published.fn \
			  bin/utils/fetch_seqs.py etc/mcra-refs.reduced.list
	$(word 1,$^) $(word 2,$^) $(word 3,$^) \
		| $(word 4,$^) $(word 5,$^) > $@

# Combine clones which have been quality trimmed with the reference sequences.
# to make the "both" file series.
# Quality trimming of the references is not required.
seq/mcra-both.luton-ampli.qtrim.fn: seq/mcra-clones.luton-ampli.qtrim.fn seq/mcra-refs.luton-ampli.fn
	cat $^ > $@

seq/mcra-both2.luton-ampli.qtrim.fn: seq/mcra-clones.luton-ampli.qtrim.fn seq/mcra-refs2.luton-ampli.fn
	cat $^ > $@

# Remove uniformative sequence {{{2
# Excise amplicon {{{3
# Search for primer hits:
res/%.psearch.out: seq/%.fn etc/primers.tsv
	primersearch -seqall $(word 1,$^) -infile $(word 2,$^) -mismatchpercent 50 -outfile $@
# I use 40% mismatchpercent so that I can include many hits before
# narrowing by other criteria.

res/%.psearch.out: seq/%.fastq etc/primers.tsv
	primersearch -seqall $(word 1,$^) -sformat fastq -infile $(word 2,$^) -mismatchpercent 50 -outfile $@
# I use 40% mismatchpercent so that I can include many hits before
# narrowing by other criteria.

res/%.psearch.tsv: bin/parse_psearch.py etc/primers.tsv res/%.psearch.out
	$(word 1,$^) $(word 2,$^) $(word 3,$^) > $@


LUTON_FIND_AMPLICON_OPTS = --primer-set luton --max-mismatch 5
# Re-orient the sequences to match the forward-reverse in etc/*.primers
# and trim to within the primers
seq/%.luton-ampli.fastq: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fastq
	$(word 1,$^) ${LUTON_FIND_AMPLICON_OPTS} --drop -f fastq -t fastq $(word 2,$^) $(word 3,$^) > $@

seq/%.luton-ampli.fn: ./bin/find_amplicon.py res/%.psearch.tsv seq/%.fn
	$^ ${LUTON_FIND_AMPLICON_OPTS} --drop > $@

# Quality trim {{{3

seq/%.qtrim.fastq: bin/qtrim_reads.py seq/%.fastq
	$^ --drop -t fastq > $@

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
seq/%.afa: seq/%.fa
	muscle < $^ > $@

seq/%.afn: bin/utils/codonalign.py seq/%.afa seq/%.frame.fn
	$^ > $@

seq/%.mcra-hmmaln.afa: etc/mcra.fungene.hmm seq/%.fa bin/utils/convert.py
	hmmalign --amino --informat FASTA $(word 1,$^) $(word 2,$^) \
		| $(word 3,$^) --in-fmt stockholm --out-fmt fasta > $@

seq/%.mcra-hmmaln.afn: bin/utils/codonalign.py seq/%.mcra-hmmaln.afa seq/%.frame.fn
	$^ > $@

seq/%.auto-refine.afa: seq/%.afa
	muscle -refine < $< > $@
# Gblocks {{{2
seq/%.gb.afn: seq/%.afn
	Gblocks $^ -t=c -p=n || [ $$? == 1 ]
	mv $^-gb $@
# gblocks always throws an error code of 1

seq/%.gb.afa: seq/%.afa
	Gblocks $^ -t=p -p=n || [ $$? == 1 ]
	mv $^-gb $@

# =======================
# Trees {{{2
tre/%.nucl.nwk: seq/%.afn
	fasttree -nt $^ > $@

tre/%.prot.nwk: seq/%.afa
	fasttree < $^ > $@

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
	pandoc -f markdown -t html5 -s --highlight-style pygments --mathjax=${MATHJAX} --toc --toc-depth=4 --css static/main.css $^ > $@

docs: ${ALL_DOCS_HTML}

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
