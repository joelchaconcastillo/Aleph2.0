# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jcc/Aleph2.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jcc/Aleph2.0/cmake

# Include any dependencies generated for this target.
include src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/flags.make

src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.o: src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/flags.make
src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.o: ../src/tools/interlevel_set_persistence_hierarchy.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.o"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.o -c /home/jcc/Aleph2.0/src/tools/interlevel_set_persistence_hierarchy.cc

src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.i"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/src/tools/interlevel_set_persistence_hierarchy.cc > CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.i

src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.s"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/src/tools/interlevel_set_persistence_hierarchy.cc -o CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.s

# Object files for target interlevel_set_persistence_hierarchy
interlevel_set_persistence_hierarchy_OBJECTS = \
"CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.o"

# External object files for target interlevel_set_persistence_hierarchy
interlevel_set_persistence_hierarchy_EXTERNAL_OBJECTS =

tools/interlevel_set_persistence_hierarchy: src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/interlevel_set_persistence_hierarchy.cc.o
tools/interlevel_set_persistence_hierarchy: src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/build.make
tools/interlevel_set_persistence_hierarchy: src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../tools/interlevel_set_persistence_hierarchy"
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/interlevel_set_persistence_hierarchy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/build: tools/interlevel_set_persistence_hierarchy

.PHONY : src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/build

src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/interlevel_set_persistence_hierarchy.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/clean

src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/src/tools /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/src/tools /home/jcc/Aleph2.0/cmake/src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/interlevel_set_persistence_hierarchy.dir/depend
