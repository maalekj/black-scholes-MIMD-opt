# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Users\Administrator\AppData\Local\Programs\CLion Nova\bin\cmake\win\x64\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Users\Administrator\AppData\Local\Programs\CLion Nova\bin\cmake\win\x64\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "D:\learning&work\mx\COE502\black-scholes-MIMD-opt"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\cmake-build-release"

# Include any dependencies generated for this target.
include CMakeFiles/black_scholes_MIMD_opt.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/black_scholes_MIMD_opt.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/black_scholes_MIMD_opt.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/black_scholes_MIMD_opt.dir/flags.make

CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj: CMakeFiles/black_scholes_MIMD_opt.dir/flags.make
CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj: D:/learning&work/mx/COE502/black-scholes-MIMD-opt/main.c
CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj: CMakeFiles/black_scholes_MIMD_opt.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir="D:\learning&work\mx\COE502\black-scholes-MIMD-opt\cmake-build-release\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj"
	C:\Users\ADMINI~1\AppData\Local\Programs\CLIONN~1\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj -MF CMakeFiles\black_scholes_MIMD_opt.dir\main.c.obj.d -o CMakeFiles\black_scholes_MIMD_opt.dir\main.c.obj -c "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\main.c"

CMakeFiles/black_scholes_MIMD_opt.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/black_scholes_MIMD_opt.dir/main.c.i"
	C:\Users\ADMINI~1\AppData\Local\Programs\CLIONN~1\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\main.c" > CMakeFiles\black_scholes_MIMD_opt.dir\main.c.i

CMakeFiles/black_scholes_MIMD_opt.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/black_scholes_MIMD_opt.dir/main.c.s"
	C:\Users\ADMINI~1\AppData\Local\Programs\CLIONN~1\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\main.c" -o CMakeFiles\black_scholes_MIMD_opt.dir\main.c.s

# Object files for target black_scholes_MIMD_opt
black_scholes_MIMD_opt_OBJECTS = \
"CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj"

# External object files for target black_scholes_MIMD_opt
black_scholes_MIMD_opt_EXTERNAL_OBJECTS =

black_scholes_MIMD_opt.exe: CMakeFiles/black_scholes_MIMD_opt.dir/main.c.obj
black_scholes_MIMD_opt.exe: CMakeFiles/black_scholes_MIMD_opt.dir/build.make
black_scholes_MIMD_opt.exe: CMakeFiles/black_scholes_MIMD_opt.dir/linkLibs.rsp
black_scholes_MIMD_opt.exe: CMakeFiles/black_scholes_MIMD_opt.dir/objects1.rsp
black_scholes_MIMD_opt.exe: CMakeFiles/black_scholes_MIMD_opt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir="D:\learning&work\mx\COE502\black-scholes-MIMD-opt\cmake-build-release\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable black_scholes_MIMD_opt.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\black_scholes_MIMD_opt.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/black_scholes_MIMD_opt.dir/build: black_scholes_MIMD_opt.exe
.PHONY : CMakeFiles/black_scholes_MIMD_opt.dir/build

CMakeFiles/black_scholes_MIMD_opt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\black_scholes_MIMD_opt.dir\cmake_clean.cmake
.PHONY : CMakeFiles/black_scholes_MIMD_opt.dir/clean

CMakeFiles/black_scholes_MIMD_opt.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "D:\learning&work\mx\COE502\black-scholes-MIMD-opt" "D:\learning&work\mx\COE502\black-scholes-MIMD-opt" "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\cmake-build-release" "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\cmake-build-release" "D:\learning&work\mx\COE502\black-scholes-MIMD-opt\cmake-build-release\CMakeFiles\black_scholes_MIMD_opt.dir\DependInfo.cmake" "--color=$(COLOR)"
.PHONY : CMakeFiles/black_scholes_MIMD_opt.dir/depend
