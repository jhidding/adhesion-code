input_files = adhesion_example.md
build_dir = ./build

compile = g++
compile_flags = -Wall -Isrc -I${HOME}/.local/include

link = g++
link_flags = -lfftw3 -lyaml-cpp

format = markdown+fenced_code_attributes
pd_call = pandoc -f $(format) --lua-filter "scripts/$(1).lua" -t plain
pd_list = $(call pd_call,list)
pd_tangle = $(call pd_call,tangle)

sources = $(shell $(pd_list) $(input_files))

cc_files = $(filter-out src/main.cc, $(filter src/%.cc, $(sources)))
obj_files = $(cc_files:%.cc=$(build_dir)/%.o)
dep_files = $(obj_files:%.o=%.d)

main_obj_file = build/src/main.o
main_dep_file = build/src/main.d

tests_cc_files = $(filter tests/%.cc, $(sources))
tests_obj_files = $(tests_cc_files:%.cc=$(build_dir)/%.o)
tests_dep_files = $(tests_obj_files:%.o=%.d)

adhesion: $(build_dir)/adhesion

run-tests: $(build_dir)/run-tests

$(sources): $(input_files)
	$(pd_tangle) $< | bash

# Include all .d files
-include $(dep_files)
-include $(tests_dep_files)
-include $(main_dep_file)

$(build_dir)/%.o : %.cc
	mkdir -p $(@D)
	$(compile) $(compile_flags) -MMD -c $< -o $@

$(build_dir)/adhesion : $(obj_files) $(main_obj_file)
	mkdir -p $(@D)
	$(link) $(link_flags) $^ -o $@

$(build_dir)/run-tests : $(obj_files) $(tests_obj_files)
	mkdir -p $(@D)
	$(link) $^ $(link_flags) -lgtest -lgmock -lpthread -o $@

.PHONY: clean
.SILENT: test

clean:
	-rm $(sources)
	-rm -r $(build_dir)

test: adhesion run-tests
	./build/run-tests
