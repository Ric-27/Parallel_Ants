# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build

# Include any dependencies generated for this target.
include CMakeFiles/spaceball.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/spaceball.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/spaceball.dir/flags.make

CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o: CMakeFiles/spaceball.dir/flags.make
CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o: /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/spaceball.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o   -c /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/spaceball.c

CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/spaceball.c > CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.i

CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/spaceball.c -o CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.s

CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.requires:

.PHONY : CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.requires

CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.provides: CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.requires
	$(MAKE) -f CMakeFiles/spaceball.dir/build.make CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.provides.build
.PHONY : CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.provides

CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.provides.build: CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o


CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o: CMakeFiles/spaceball.dir/flags.make
CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o: /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/vmath.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o   -c /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/vmath.c

CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/vmath.c > CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.i

CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0/progs/demos/spaceball/vmath.c -o CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.s

CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.requires:

.PHONY : CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.requires

CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.provides: CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.requires
	$(MAKE) -f CMakeFiles/spaceball.dir/build.make CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.provides.build
.PHONY : CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.provides

CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.provides.build: CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o


# Object files for target spaceball
spaceball_OBJECTS = \
"CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o" \
"CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o"

# External object files for target spaceball
spaceball_EXTERNAL_OBJECTS =

bin/spaceball: CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o
bin/spaceball: CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o
bin/spaceball: CMakeFiles/spaceball.dir/build.make
bin/spaceball: /usr/lib/x86_64-linux-gnu/libGLU.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libGL.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libSM.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libICE.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libX11.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXext.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXrandr.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXi.so
bin/spaceball: lib/libglut.so.3.10.0
bin/spaceball: /usr/lib/x86_64-linux-gnu/libGL.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libSM.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libICE.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libX11.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXext.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXrandr.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
bin/spaceball: /usr/lib/x86_64-linux-gnu/libXi.so
bin/spaceball: CMakeFiles/spaceball.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable bin/spaceball"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/spaceball.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/spaceball.dir/build: bin/spaceball

.PHONY : CMakeFiles/spaceball.dir/build

CMakeFiles/spaceball.dir/requires: CMakeFiles/spaceball.dir/progs/demos/spaceball/spaceball.c.o.requires
CMakeFiles/spaceball.dir/requires: CMakeFiles/spaceball.dir/progs/demos/spaceball/vmath.c.o.requires

.PHONY : CMakeFiles/spaceball.dir/requires

CMakeFiles/spaceball.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/spaceball.dir/cmake_clean.cmake
.PHONY : CMakeFiles/spaceball.dir/clean

CMakeFiles/spaceball.dir/depend:
	cd /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0 /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/freeglut-3.0.0 /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build /home/ric/Documents/ENSTA/2/IN203/Trabajo_Propio/Projet_laby_ants/libgui/thirdparty/build/CMakeFiles/spaceball.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/spaceball.dir/depend

