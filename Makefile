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

download-seqs:
	cd raw ; \
	curl -u $$SEQCORE_USER:$$SEQCORE_PSWD -O "http://seqcore.brcf.med.umich.edu/users/schmidt/smith/[2476527-2476582].M13REV.ab1"
	rm -f raw/2476550.M13REV.seq raw/2476529.M13REV.seq
# These two don't exist on the server.
# How do I remove them automatically?

download-traces:
	cd raw/traces ; \
	curl -u $$SEQCORE_USER:$$SEQCORE_PSWD -O "http://seqcore.brcf.med.umich.edu/users/schmidt/smith/[2476527-2476582].M13REV.ab1"
	rm -f raw/2476550.M13REV.ab1 raw/2476529.M13REV.ab1

seq/clones.raw-names.with-suspect.fn:
	echo "" > $@
	for seq_file in raw/*.seq; do \
		base=$$(basename $$seq_file) ; \
		echo ">$${base/.seq/}" >> $@ ; \
		cat $$seq_file >> $@ ; \
	done
	sed -i '1,1d' $@

seq/clones.with-suspect.fn: bin/utils/rename_seqs.py etc/clones.names.tsv seq/clones.raw-names.with-suspect.fn
	$^ > $@

seq/clones.fn: bin/utils/drop_seqs.py etc/clones.suspect.list seq/clones.with-suspect.fn
	$^ > $@

seq/refs.fn: bin/utils/rename_seqs.py etc/refs.names.tsv raw/mcra.refs.fn
	$^ > $@

seq/clones.with-refs.fn: seq/clones.fn seq/refs.fn
	cat $^ > $@

res/%.primersearch-luton.out: seq/%.fn etc/luton_primers.txt
	primersearch \
		-seqall seq/$*.fn \
		-infile etc/luton_primers.txt \
		-mismatchpercent 40 \
		-outfile $@
# I use 40% mismatchpercent so that I can include many hits before
# narrowing by other criteria.

res/%.primersearch-luton.tsv: res/%.primersearch-luton.out bin/mung_primersearch.py
	bin/mung_primersearch.py < $< > $@

seq/%.qc.fn: bin/clean_clones.py res/%.primersearch-luton.tsv seq/%.fn
	$^ GGTGG >$@
# This step is certainly not universal, and causes me to lose two of my
# references.
# I think that I lose too many sequences on this step.  Will a Quality Score
# based cleanup be better?

# Translation:
seq/%.fa: bin/utils/translate.py seq/%.fn
	cat seq/$*.fn \
		| bin/utils/translate.py \
		> $@

# Alignment
seq/%.afa: seq/%.fa
	cat $^ \
		| muscle \
		> $@

# Reordering:  *`muscle` puts sequences out of order
seq/%.ord.afa: bin/utils/ls_ids.py bin/utils/fetch_seqs.py seq/%.afa seq/%.fn
	$(eval $*_TMP := $(shell mktemp))
	bin/utils/ls_ids.py seq/$*.fn > ${$*_TMP}
	bin/utils/fetch_seqs.py ${$*_TMP} seq/$*.afa > $@
	rm ${$*_TMP}


# Backalign:
seq/%.afn: bin/utils/backalign.py seq/%.ord.afa seq/%.fn
	$^ > $@


# Trees
tre/%.nucl.nwk: seq/%.afn
	fasttree -nt $^ > $@

tre/%.prot.nwk: seq/%.afa
	fasttree < $^ > $@


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
#  Documentation Recipes {{{1
# =======================
ALL_DOCS = TEMPLATE NOTE
ALL_DOCS_HTML = $(addsuffix .html,${ALL_DOCS})
MATHJAX = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

%.html: %.md
	cat $^ \
	| pandoc -f markdown -t html5 -s\
				--highlight-style pygments \
				--mathjax=${MATHJAX} \
				--toc --toc-depth=4 \
				--css static/main.css \
	> $@

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
