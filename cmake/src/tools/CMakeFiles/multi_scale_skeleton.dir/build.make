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
include src/tools/CMakeFiles/multi_scale_skeleton.dir/depend.make

# Include the progress variables for this target.
include src/tools/CMakeFiles/multi_scale_skeleton.dir/progress.make

# Include the compile flags for this target's objects.
include src/tools/CMakeFiles/multi_scale_skeleton.dir/flags.make

src/tools/CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.o: src/tools/CMakeFiles/multi_scale_skeleton.dir/flags.make
src/tools/CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.o: ../src/tools/multi_scale_skeleton.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/tools/CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.o"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.o -c /home/jcc/Aleph2.0/src/tools/multi_scale_skeleton.cc

src/tools/CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.i"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jcc/Aleph2.0/src/tools/multi_scale_skeleton.cc > CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.i

src/tools/CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.s"
	cd /home/jcc/Aleph2.0/cmake/src/tools && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jcc/Aleph2.0/src/tools/multi_scale_skeleton.cc -o CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.s

# Object files for target multi_scale_skeleton
multi_scale_skeleton_OBJECTS = \
"CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.o"

# External object files for target multi_scale_skeleton
multi_scale_skeleton_EXTERNAL_OBJECTS =

tools/multi_scale_skeleton: src/tools/CMakeFiles/multi_scale_skeleton.dir/multi_scale_skeleton.cc.o
tools/multi_scale_skeleton: src/tools/CMakeFiles/multi_scale_skeleton.dir/build.make
tools/multi_scale_skeleton: src/tools/CMakeFiles/multi_scale_skeleton.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jcc/Aleph2.0/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../tools/multi_scale_skeleton"
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/multi_scale_skeleton.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/tools/CMakeFiles/multi_scale_skeleton.dir/build: tools/multi_scale_skeleton

.PHONY : src/tools/CMakeFiles/multi_scale_skeleton.dir/build

src/tools/CMakeFiles/multi_scale_skeleton.dir/clean:
	cd /home/jcc/Aleph2.0/cmake/src/tools && $(CMAKE_COMMAND) -P CMakeFiles/multi_scale_skeleton.dir/cmake_clean.cmake
.PHONY : src/tools/CMakeFiles/multi_scale_skeleton.dir/clean

src/tools/CMakeFiles/multi_scale_skeleton.dir/depend:
	cd /home/jcc/Aleph2.0/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jcc/Aleph2.0 /home/jcc/Aleph2.0/src/tools /home/jcc/Aleph2.0/cmake /home/jcc/Aleph2.0/cmake/src/tools /home/jcc/Aleph2.0/cmake/src/tools/CMakeFiles/multi_scale_skeleton.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/tools/CMakeFiles/multi_scale_skeleton.dir/depend
