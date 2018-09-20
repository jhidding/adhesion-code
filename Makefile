input_files = adhesion_example.md
build_dir = ./build

hdf5_cflags = $(shell pkg-config --cflags hdf5)
hdf5_libs = $(shell pkg-config --libs hdf5) -lhdf5_cpp

cgal_cflags = -frounding-math
cgal_libs = -lm -lCGAL -lgmp -lboost_thread -lmpfr

gsl_libs = $(shell pkg-config --libs gsl)

coverage_cflags = --coverage

compile = g++
compile_flags = -std=c++17 -O3 -Wall -Isrc -I${HOME}/.local/include $(hdf5_cflags) $(cgal_cflags)

link = g++
link_flags = -lfftw3 -lyaml-cpp -lfmt $(hdf5_libs) $(cgal_libs) $(gsl_libs)

# ===========================================================================

SHELL := /bin/bash

format = markdown+fenced_code_attributes+citations+all_symbols_escapable
pandoc_filters = pandoc-eqnos pandoc-fignos pandoc-citeproc
report_args = --toc $(pandoc_filters:%=--filter %) --lua-filter "scripts/annotate-code-blocks.lua" --template scripts/eisvogel.tex --listings 
html_args = -s --toc --toc-depth=2 $(pandoc_filters:%=--filter %) --lua-filter "scripts/annotate-code-blocks.lua" --mathjax --css "style.css"

pd_call = pandoc -f $(format) --lua-filter "scripts/$(1).lua" -t plain
pd_list = $(call pd_call,list)
pd_tangle = $(call pd_call,tangle)

pdf_files = $(input_files:%.md=$(build_dir)/%.pdf)
html_files = $(input_files:%.md=$(build_dir)/html/%.html)
sources = $(shell $(pd_list) $(input_files))

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

run-tests: $(build_dir)/run-tests

tangle: $(input_files)
	mkdir -p $(build_dir)
	$(pd_tangle) $< > $(build_dir)/tangle.sh
	source $(build_dir)/tangle.sh

report: $(pdf_files)

html: $(html_files)
	cp -r figures build/html
	cp scripts/style.css build/html

$(build_dir)/html/%.html: %.md
	mkdir -p $(build_dir)/html
	pandoc $< -f $(format) $(html_args) -t html5 -o $@

$(build_dir)/%.pdf : %.md
	pandoc $< -f $(format) $(report_args) -t latex -o $@ --pdf-engine=xelatex

$(sources): tangle

# Include all .d files
-include $(dep_files)
-include $(tests_dep_files)
-include $(main_dep_file)

$(build_dir)/%.o : %.cc
	mkdir -p $(@D)
	$(compile) $(compile_flags) -MMD -c $< -o $@

# Link main executable
$(build_dir)/adhesion : $(obj_files) $(main_obj_file)
	mkdir -p $(@D)
	$(link) $^ $(link_flags) -o $@

# Link testing exectuable
$(build_dir)/run-tests : $(obj_files) $(tests_obj_files)
	mkdir -p $(@D)
	$(link) $^ $(link_flags) -lgtest -lgmock -lpthread -o $@

clean:
	-rm $(sources)
	-rm -r $(build_dir)

test: adhesion run-tests
	./build/run-tests
