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
CMAKE_SOURCE_DIR = /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build

# Include any dependencies generated for this target.
include CMakeFiles/project.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/project.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/project.dir/flags.make

CMakeFiles/project.dir/src/estimatePSD.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/estimatePSD.cpp.o: ../src/estimatePSD.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/project.dir/src/estimatePSD.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/estimatePSD.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/estimatePSD.cpp

CMakeFiles/project.dir/src/estimatePSD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/estimatePSD.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/estimatePSD.cpp > CMakeFiles/project.dir/src/estimatePSD.cpp.i

CMakeFiles/project.dir/src/estimatePSD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/estimatePSD.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/estimatePSD.cpp -o CMakeFiles/project.dir/src/estimatePSD.cpp.s

CMakeFiles/project.dir/src/filter.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/filter.cpp.o: ../src/filter.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/project.dir/src/filter.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/filter.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/filter.cpp

CMakeFiles/project.dir/src/filter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/filter.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/filter.cpp > CMakeFiles/project.dir/src/filter.cpp.i

CMakeFiles/project.dir/src/filter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/filter.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/filter.cpp -o CMakeFiles/project.dir/src/filter.cpp.s

CMakeFiles/project.dir/src/fourier.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/fourier.cpp.o: ../src/fourier.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/project.dir/src/fourier.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/fourier.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/fourier.cpp

CMakeFiles/project.dir/src/fourier.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/fourier.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/fourier.cpp > CMakeFiles/project.dir/src/fourier.cpp.i

CMakeFiles/project.dir/src/fourier.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/fourier.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/fourier.cpp -o CMakeFiles/project.dir/src/fourier.cpp.s

CMakeFiles/project.dir/src/genfunc.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/genfunc.cpp.o: ../src/genfunc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/project.dir/src/genfunc.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/genfunc.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/genfunc.cpp

CMakeFiles/project.dir/src/genfunc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/genfunc.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/genfunc.cpp > CMakeFiles/project.dir/src/genfunc.cpp.i

CMakeFiles/project.dir/src/genfunc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/genfunc.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/genfunc.cpp -o CMakeFiles/project.dir/src/genfunc.cpp.s

CMakeFiles/project.dir/src/iofunc.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/iofunc.cpp.o: ../src/iofunc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/project.dir/src/iofunc.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/iofunc.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/iofunc.cpp

CMakeFiles/project.dir/src/iofunc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/iofunc.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/iofunc.cpp > CMakeFiles/project.dir/src/iofunc.cpp.i

CMakeFiles/project.dir/src/iofunc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/iofunc.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/iofunc.cpp -o CMakeFiles/project.dir/src/iofunc.cpp.s

CMakeFiles/project.dir/src/logfunc.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/logfunc.cpp.o: ../src/logfunc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/project.dir/src/logfunc.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/logfunc.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/logfunc.cpp

CMakeFiles/project.dir/src/logfunc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/logfunc.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/logfunc.cpp > CMakeFiles/project.dir/src/logfunc.cpp.i

CMakeFiles/project.dir/src/logfunc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/logfunc.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/logfunc.cpp -o CMakeFiles/project.dir/src/logfunc.cpp.s

CMakeFiles/project.dir/src/project.cpp.o: CMakeFiles/project.dir/flags.make
CMakeFiles/project.dir/src/project.cpp.o: ../src/project.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/project.dir/src/project.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project.dir/src/project.cpp.o -c /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/project.cpp

CMakeFiles/project.dir/src/project.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project.dir/src/project.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/project.cpp > CMakeFiles/project.dir/src/project.cpp.i

CMakeFiles/project.dir/src/project.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project.dir/src/project.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/src/project.cpp -o CMakeFiles/project.dir/src/project.cpp.s

# Object files for target project
project_OBJECTS = \
"CMakeFiles/project.dir/src/estimatePSD.cpp.o" \
"CMakeFiles/project.dir/src/filter.cpp.o" \
"CMakeFiles/project.dir/src/fourier.cpp.o" \
"CMakeFiles/project.dir/src/genfunc.cpp.o" \
"CMakeFiles/project.dir/src/iofunc.cpp.o" \
"CMakeFiles/project.dir/src/logfunc.cpp.o" \
"CMakeFiles/project.dir/src/project.cpp.o"

# External object files for target project
project_EXTERNAL_OBJECTS =

project: CMakeFiles/project.dir/src/estimatePSD.cpp.o
project: CMakeFiles/project.dir/src/filter.cpp.o
project: CMakeFiles/project.dir/src/fourier.cpp.o
project: CMakeFiles/project.dir/src/genfunc.cpp.o
project: CMakeFiles/project.dir/src/iofunc.cpp.o
project: CMakeFiles/project.dir/src/logfunc.cpp.o
project: CMakeFiles/project.dir/src/project.cpp.o
project: CMakeFiles/project.dir/build.make
project: CMakeFiles/project.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable project"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/project.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/project.dir/build: project

.PHONY : CMakeFiles/project.dir/build

CMakeFiles/project.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/project.dir/cmake_clean.cmake
.PHONY : CMakeFiles/project.dir/clean

CMakeFiles/project.dir/depend:
	cd /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build /home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/CMakeFiles/project.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/project.dir/depend
