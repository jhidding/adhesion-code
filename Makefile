# Input files
input_files = adhesion_example.md python_bindings.md appendix.md

# Target directories
build_dir = ./build
html_dir = ./docs

# Target files
pdf_files = $(build_dir)/main.pdf
html_files = $(html_dir)/index.html

static_files := $(shell find -L figures -type f)
static_targets := $(static_files:%=docs/%)

# Get HDF5 compilation flags using pkg-config
# hdf5_cflags = $(shell pkg-config --cflags hdf5)
# hdf5_libs = $(shell pkg-config --libs hdf5) -lhdf5_cpp

# Get HDF5 compilation flags manually if pkg-config does not work
hdf5_cc := $(shell which h5c++)
hdf5_root := $(shell dirname `dirname $(hdf5_cc)`)
hdf5_cflags := $(shell h5c++ -show '%' | cut -d% -f1 | cut -d' ' -f2-) -I$(hdf5_root)/include
hdf5_libs := $(shell h5c++ -show '%' | cut -d% -f2 | cut -d' ' -f2-)

# CGAL compile flags when compiling with GCC
# Note that CLANG does not support -frounding-math
cgal_cflags = -frounding-math
cgal_libs = -lm -lgmp -lboost_thread -lmpfr

# Uncomment if you want to compile with TBB
# tbb_cflags = -DCGAL_LINKED_WITH_TBB -pthread
# tbb_libs = -ltbb -latomic -ltbbmalloc -pthread

gsl_libs := $(shell pkg-config --libs gsl)

coverage_cflags = --coverage

# In case we're compiling with GCC 6
# cflags = -D_GLIBCXX_USE_CXX11_ABI=0 -I${HOME}/.local/include $(tbb_cflags)
# cflags = -O3
cflags = -O3 -Iinclude -I${HOME}/.local/include

# If some of the dependencies are installed locally
libs = -L${HOME}/.local/lib
# libs =

sources := $(shell entangled list)

cc_files := $(filter-out src/main.cc, $(filter src/%.cc, $(sources)))
obj_files = $(cc_files:%.cc=$(build_dir)/%.o)
dep_files = $(obj_files:%.o=%.d)

python_sources = $(cc_files)
python_build_dir = $(build_dir)/python-module
python_cflags := -fPIC $(shell python -m pybind11 --includes)
# python_libs = $(shell pkg-config --libs python3)
python_obj_files = $(python_sources:%.cc=$(python_build_dir)/%.o)
python_dep_files = $(python_sources:%.cc=$(python_build_dir)/%.d)
python_target := $(python_build_dir)/adhesion$(shell python3-config --extension-suffix)

compile = g++
compile_flags = -std=c++17 -Wall -Isrc  $(cflags) $(hdf5_cflags) $(cgal_cflags) $(tbb_cflags)
link = g++
link_flags = -lstdc++fs -lfftw3 -lyaml-cpp -lfmt $(libs) $(hdf5_libs) $(cgal_libs) $(gsl_libs) $(tbb_libs)

# ===========================================================================

SHELL := /bin/bash

format = markdown+fenced_code_attributes+citations+all_symbols_escapable+fenced_divs+multiline_tables
pandoc_filters = pandoc-eqnos pandoc-fignos pandoc-citeproc
pandoc_latex_filters = pandoc-eqnos pandoc-fignos
report_args = --toc $(pandoc_filters:%=--filter %) --lua-filter "scripts/annotate-code-blocks.lua" --template scripts/eisvogel.tex --listings
latex_args = --toc $(pandoc_latex_filters:%=--filter %) --lua-filter "scripts/annotate-code-blocks.lua" --natbib --listings

main_obj_file = build/src/main.o
main_dep_file = build/src/main.d

tests_cc_files = $(filter tests/%.cc, $(sources))
tests_obj_files = $(tests_cc_files:%.cc=$(build_dir)/%.o)
tests_dep_files = $(tests_obj_files:%.o=%.d)

# Arguments to Pandoc; these are reasonable defaults
pandoc_args += --template bootstrap/template.html
pandoc_args += --css css/mods.css
pandoc_args += -t html5 -s --mathjax --toc
pandoc_args += --toc-depth 2
pandoc_args += --filter pandoc-bootstrap
pandoc_args += --filter pandoc-eqnos
pandoc_args += --filter pandoc-fignos
pandoc_args += --filter pandoc-citeproc
pandoc_args += -f markdown+multiline_tables+simple_tables
pandoc_args += --highlight-style tango

# ===========================================================================

.PHONY: clean tangle adhesion run-tests site python-module
# .SILENT: test

watch:
	@tmux new-session make --no-print-directory watch-pandoc \; \
		split-window -v make --no-print-directory watch-browser-sync \; \
		split-window -v entangled daemon \; \
		select-layout even-vertical \;

watch-pandoc:
	@while true; do \
		inotifywait -e close_write bootstrap Makefile $(input_files); \
		make site; \
	done

watch-browser-sync:
	browser-sync start -w -s docs

site: docs/index.html docs/css/mods-scrollspy.css $(static_targets)

docs/index.html: $(input_files) Makefile bootstrap/template-scrollspy.html
	@mkdir -p docs
	pandoc $(pandoc_args) $(input_files) -o $@

docs/css/mods-scrollspy.css: bootstrap/mods-scrollspy.css
	@mkdir -p docs/css
	cp $< $@

$(static_targets): docs/%: %
	@mkdir -p $(dir $@)
	cp $< $@

# ==================================================================

adhesion: $(build_dir)/adhesion

python-module: $(python_target)

parallel-test: $(build_dir)/parallel-test

run-tests: $(build_dir)/run-tests

report: $(pdf_files)

latex: $(build_dir)/adhesion_example.tex
	cat $(build_dir)/adhesion_example.tex \
		| perl -0777 -pe 's/(\\emph\{.*\})\n\n\\begin\{lstlisting\}/\\par\\needspace{4\\baselineskip}\\vspace{8pt}$$1\\vspace{-8pt}\\begin{lstlisting}/g' \
		| perl -0777 -pe 's/\\includegraphics\{/\\includegraphics[width=\\textwidth]{/g' \
		| perl -0777 -pe 's/\\includegraphics\[(.*)\]\{(.*)\.svg\}/\\includegraphics[$$1]{$$2}/g' \
		| perl -0777 -pe 's/\\includegraphics\[(.*),height=.*\]\{(.*)\}/\\includegraphics[$$1]{$$2}/g' \
		> $(build_dir)/adhesion_code.tex

$(build_dir)/%.tex : %.md
	pandoc $^ -f $(format) $(latex_args) -t latex -o $@

$(build_dir)/%.pdf : $(input_files)
	pandoc $^ -f $(format) $(report_args) -t latex -o $@ --pdf-engine=xelatex

# Include all .d files
-include $(dep_files)
-include $(tests_dep_files)
-include $(main_dep_file)
-include $(python_dep_files)

$(python_build_dir)/%.o: %.cc Makefile
	@echo [Compiling]: $<
	@mkdir -p $(@D)
	$(compile) $(compile_flags) $(python_cflags) -fPIC -MMD -c $< -o $@

$(python_target): $(python_obj_files)
	$(link) $^ $(link_flags) -shared -o $@

$(build_dir)/%.o : %.cc Makefile
	@echo [Compiling]: $<
	@mkdir -p $(@D)
	@$(compile) $(compile_flags) -MMD -c $< -o $@

$(build_dir)/parallel-test : $(build_dir)/examples/cgal-parallel.o
	@mkdir -p $(@D)
	$(link) $^ $(link_flags) -o $@

# Link main executable
$(build_dir)/adhesion : $(obj_files) $(main_obj_file)
	@echo [Linking] $@: $^
	@mkdir -p $(@D)
	$(link) $^ $(link_flags) -o $@

# Link testing exectuable
$(build_dir)/run-tests : $(obj_files) $(tests_obj_files)
	@mkdir -p $(@D)
	$(link) $^ $(link_flags) -lgtest -lgmock -lpthread -o $@

clean:
	-rm -r $(build_dir)

test: adhesion run-tests
	./build/run-tests
