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
include examples/CMakeFiles/image_analysis.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/image_analysis.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/image_analysis.dir/flags.make

examples/CMakeFiles/image_analysis.dir/image_analysis.cc.o: examples/CMakeFiles/image_analysis.dir/flags.make
examples/CMakeFiles/image_analysis.dir/image_analysis.cc.o: ../examples/image_analysis.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/image_analysis.dir/image_analysis.cc.o"
	cd /home/jcc/Aleph2.0/cmake/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/image_analysis.dir/image_analysis.cc.o -c /home/jcc/Aleph2.0/examples/image_analysis.cc

examples/CMakeFiles/image_analysis.dir/image_analysis.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/image_analysis.dir/image_analysis.cc.i"
	cd /home/jcc/Aleph2.0/cmake/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/examples/image_analysis.cc > CMakeFiles/image_analysis.dir/image_analysis.cc.i

examples/CMakeFiles/image_analysis.dir/image_analysis.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/image_analysis.dir/image_analysis.cc.s"
	cd /home/jcc/Aleph2.0/cmake/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/examples/image_analysis.cc -o CMakeFiles/image_analysis.dir/image_analysis.cc.s

# Object files for target image_analysis
image_analysis_OBJECTS = \
"CMakeFiles/image_analysis.dir/image_analysis.cc.o"

# External object files for target image_analysis
image_analysis_EXTERNAL_OBJECTS =

examples/image_analysis: examples/CMakeFiles/image_analysis.dir/image_analysis.cc.o
examples/image_analysis: examples/CMakeFiles/image_analysis.dir/build.make
examples/image_analysis: examples/CMakeFiles/image_analysis.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable image_analysis"
	cd /home/jcc/Aleph2.0/cmake/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/image_analysis.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/image_analysis.dir/build: examples/image_analysis

.PHONY : examples/CMakeFiles/image_analysis.dir/build

examples/CMakeFiles/image_analysis.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/examples && $(CMAKE_COMMAND) -P CMakeFiles/image_analysis.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/image_analysis.dir/clean

examples/CMakeFiles/image_analysis.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/examples /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/examples /home/jcc/Aleph2.0/cmake/examples/CMakeFiles/image_analysis.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/image_analysis.dir/depend
