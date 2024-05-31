# Input files
input_files = adhesion_example.md appendix.md

# Target directories
build_dir = ./build
html_dir = ./docs

# Target files
pdf_files = $(build_dir)/main.pdf
html_files = $(html_dir)/index.html

# Get HDF5 compilation flags using pkg-config
hdf5_cflags = $(shell pkg-config --cflags hdf5)
hdf5_libs = $(shell pkg-config --libs hdf5) -lhdf5_cpp

# Get HDF5 compilation flags manually if pkg-config does not work
# hdf5_cc = $(shell which h5c++)
# hdf5_root = $(shell dirname `dirname $(hdf5_cc)`)
# hdf5_cflags = $(shell h5c++ -show '%' | cut -d% -f1 | cut -d' ' -f2-) -I$(hdf5_root)/include
# hdf5_libs = $(shell h5c++ -show '%' | cut -d% -f2 | cut -d' ' -f2-)

# CGAL compile flags when compiling with GCC
# Note that CLANG does not support -frounding-math
cgal_cflags = -frounding-math
cgal_libs = -lm -lgmp -lboost_thread -lmpfr

# Uncomment if you want to compile with TBB
# tbb_cflags = -DCGAL_LINKED_WITH_TBB -pthread
# tbb_libs = -ltbb -latomic -ltbbmalloc -pthread

gsl_libs = $(shell pkg-config --libs gsl)

coverage_cflags = --coverage

# In case we're compiling with GCC 6
# cflags = -D_GLIBCXX_USE_CXX11_ABI=0 -I${HOME}/.local/include $(tbb_cflags)
# cflags = -O3
cflags = -O3 -g -Iinclude 

# If some of the dependencies are installed locally
# libs = -L${HOME}/.local/lib
# libs =

compile = g++
# compile = clang
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
html_args = -s --toc --toc-depth=3 $(pandoc_filters:%=--filter %) --lua-filter "scripts/annotate-code-blocks.lua" --mathjax --css "style.css" --base-header-level=2

sources = $(shell find src -name '*.cc')

cc_files = $(filter-out src/main.cc, $(filter src/%.cc, $(sources)))
obj_files = $(cc_files:%.cc=$(build_dir)/%.o)
dep_files = $(obj_files:%.o=%.d)

main_obj_file = build/src/main.o
main_dep_file = build/src/main.d

tests_cc_files = $(filter tests/%.cc, $(sources))
tests_obj_files = $(tests_cc_files:%.cc=$(build_dir)/%.o)
tests_dep_files = $(tests_obj_files:%.o=%.d)

# ===========================================================================

.PHONY: clean tangle adhesion run-tests
# .SILENT: test

adhesion: $(build_dir)/adhesion

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

html: $(html_files)
	cp -r figures $(html_dir)
	cp scripts/style.css $(html_dir)

$(html_dir)/%.html: $(input_files)
	mkdir -p $(html_dir)
	pandoc $^ -f $(format) $(html_args) -t html5 -o $@

$(build_dir)/%.pdf : $(input_files)
	pandoc $^ -f $(format) $(report_args) -t latex -o $@ --pdf-engine=xelatex

# Include all .d files
-include $(dep_files)
-include $(tests_dep_files)
-include $(main_dep_file)

$(build_dir)/%.o : %.cc Makefile
	@mkdir -p $(@D)
	$(compile) $(compile_flags) -MMD -c $< -o $@

$(build_dir)/parallel-test : $(build_dir)/examples/cgal-parallel.o
	@mkdir -p $(@D)
	$(link) $^ $(link_flags) -o $@

# Link main executable
$(build_dir)/adhesion : $(obj_files) $(main_obj_file)
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
