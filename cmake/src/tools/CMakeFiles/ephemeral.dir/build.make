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
include src/tools/CMakeFiles/ephemeral.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/ephemeral.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/ephemeral.dir/flags.make

src/tools/CMakeFiles/ephemeral.dir/ephemeral.cc.o: src/tools/CMakeFiles/ephemeral.dir/flags.make
src/tools/CMakeFiles/ephemeral.dir/ephemeral.cc.o: ../src/tools/ephemeral.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/tools/CMakeFiles/ephemeral.dir/ephemeral.cc.o"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ephemeral.dir/ephemeral.cc.o -c /home/jcc/Aleph2.0/src/tools/ephemeral.cc

src/tools/CMakeFiles/ephemeral.dir/ephemeral.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ephemeral.dir/ephemeral.cc.i"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/src/tools/ephemeral.cc > CMakeFiles/ephemeral.dir/ephemeral.cc.i

src/tools/CMakeFiles/ephemeral.dir/ephemeral.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ephemeral.dir/ephemeral.cc.s"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/src/tools/ephemeral.cc -o CMakeFiles/ephemeral.dir/ephemeral.cc.s

# Object files for target ephemeral
ephemeral_OBJECTS = \
"CMakeFiles/ephemeral.dir/ephemeral.cc.o"

# External object files for target ephemeral
ephemeral_EXTERNAL_OBJECTS =

tools/ephemeral: src/tools/CMakeFiles/ephemeral.dir/ephemeral.cc.o
tools/ephemeral: src/tools/CMakeFiles/ephemeral.dir/build.make
tools/ephemeral: src/tools/CMakeFiles/ephemeral.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../tools/ephemeral"
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ephemeral.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/ephemeral.dir/build: tools/ephemeral

.PHONY : src/tools/CMakeFiles/ephemeral.dir/build

src/tools/CMakeFiles/ephemeral.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/ephemeral.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/ephemeral.dir/clean

src/tools/CMakeFiles/ephemeral.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/src/tools /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/src/tools /home/jcc/Aleph2.0/cmake/src/tools/CMakeFiles/ephemeral.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/ephemeral.dir/depend

