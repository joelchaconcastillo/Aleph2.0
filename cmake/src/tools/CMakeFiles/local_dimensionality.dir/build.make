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
include src/tools/CMakeFiles/local_dimensionality.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/local_dimensionality.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/local_dimensionality.dir/flags.make

src/tools/CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.o: src/tools/CMakeFiles/local_dimensionality.dir/flags.make
src/tools/CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.o: ../src/tools/local_dimensionality.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/tools/CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.o"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.o -c /home/jcc/Aleph2.0/src/tools/local_dimensionality.cc

src/tools/CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.i"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/src/tools/local_dimensionality.cc > CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.i

src/tools/CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.s"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/src/tools/local_dimensionality.cc -o CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.s

# Object files for target local_dimensionality
local_dimensionality_OBJECTS = \
"CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.o"

# External object files for target local_dimensionality
local_dimensionality_EXTERNAL_OBJECTS =

tools/local_dimensionality: src/tools/CMakeFiles/local_dimensionality.dir/local_dimensionality.cc.o
tools/local_dimensionality: src/tools/CMakeFiles/local_dimensionality.dir/build.make
tools/local_dimensionality: src/tools/CMakeFiles/local_dimensionality.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../tools/local_dimensionality"
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/local_dimensionality.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/local_dimensionality.dir/build: tools/local_dimensionality

.PHONY : src/tools/CMakeFiles/local_dimensionality.dir/build

src/tools/CMakeFiles/local_dimensionality.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/local_dimensionality.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/local_dimensionality.dir/clean

src/tools/CMakeFiles/local_dimensionality.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/src/tools /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/src/tools /home/jcc/Aleph2.0/cmake/src/tools/CMakeFiles/local_dimensionality.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/local_dimensionality.dir/depend
