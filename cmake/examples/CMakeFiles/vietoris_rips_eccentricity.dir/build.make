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
include examples/CMakeFiles/vietoris_rips_eccentricity.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/vietoris_rips_eccentricity.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/vietoris_rips_eccentricity.dir/flags.make

examples/CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.o: examples/CMakeFiles/vietoris_rips_eccentricity.dir/flags.make
examples/CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.o: ../examples/vietoris_rips_eccentricity.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.o"
	cd /home/jcc/Aleph2.0/cmake/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.o -c /home/jcc/Aleph2.0/examples/vietoris_rips_eccentricity.cc

examples/CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.i"
	cd /home/jcc/Aleph2.0/cmake/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/examples/vietoris_rips_eccentricity.cc > CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.i

examples/CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.s"
	cd /home/jcc/Aleph2.0/cmake/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/examples/vietoris_rips_eccentricity.cc -o CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.s

# Object files for target vietoris_rips_eccentricity
vietoris_rips_eccentricity_OBJECTS = \
"CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.o"

# External object files for target vietoris_rips_eccentricity
vietoris_rips_eccentricity_EXTERNAL_OBJECTS =

examples/vietoris_rips_eccentricity: examples/CMakeFiles/vietoris_rips_eccentricity.dir/vietoris_rips_eccentricity.cc.o
examples/vietoris_rips_eccentricity: examples/CMakeFiles/vietoris_rips_eccentricity.dir/build.make
examples/vietoris_rips_eccentricity: examples/CMakeFiles/vietoris_rips_eccentricity.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vietoris_rips_eccentricity"
	cd /home/jcc/Aleph2.0/cmake/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vietoris_rips_eccentricity.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/vietoris_rips_eccentricity.dir/build: examples/vietoris_rips_eccentricity

.PHONY : examples/CMakeFiles/vietoris_rips_eccentricity.dir/build

examples/CMakeFiles/vietoris_rips_eccentricity.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/examples && $(CMAKE_COMMAND) -P CMakeFiles/vietoris_rips_eccentricity.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/vietoris_rips_eccentricity.dir/clean

examples/CMakeFiles/vietoris_rips_eccentricity.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/examples /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/examples /home/jcc/Aleph2.0/cmake/examples/CMakeFiles/vietoris_rips_eccentricity.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/vietoris_rips_eccentricity.dir/depend
