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
include src/tools/CMakeFiles/persistence_diagram_statistics.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/persistence_diagram_statistics.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/persistence_diagram_statistics.dir/flags.make

src/tools/CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.o: src/tools/CMakeFiles/persistence_diagram_statistics.dir/flags.make
src/tools/CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.o: ../src/tools/persistence_diagram_statistics.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/tools/CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.o"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.o -c /home/jcc/Aleph2.0/src/tools/persistence_diagram_statistics.cc

src/tools/CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.i"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/src/tools/persistence_diagram_statistics.cc > CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.i

src/tools/CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.s"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/src/tools/persistence_diagram_statistics.cc -o CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.s

# Object files for target persistence_diagram_statistics
persistence_diagram_statistics_OBJECTS = \
"CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.o"

# External object files for target persistence_diagram_statistics
persistence_diagram_statistics_EXTERNAL_OBJECTS =

tools/persistence_diagram_statistics: src/tools/CMakeFiles/persistence_diagram_statistics.dir/persistence_diagram_statistics.cc.o
tools/persistence_diagram_statistics: src/tools/CMakeFiles/persistence_diagram_statistics.dir/build.make
tools/persistence_diagram_statistics: src/tools/CMakeFiles/persistence_diagram_statistics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../tools/persistence_diagram_statistics"
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/persistence_diagram_statistics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/persistence_diagram_statistics.dir/build: tools/persistence_diagram_statistics

.PHONY : src/tools/CMakeFiles/persistence_diagram_statistics.dir/build

src/tools/CMakeFiles/persistence_diagram_statistics.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/persistence_diagram_statistics.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/persistence_diagram_statistics.dir/clean

src/tools/CMakeFiles/persistence_diagram_statistics.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/src/tools /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/src/tools /home/jcc/Aleph2.0/cmake/src/tools/CMakeFiles/persistence_diagram_statistics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/persistence_diagram_statistics.dir/depend

