# I think I understand why meta-build systems exist now. Guess I should start
# using CMake or Premake or something.

# Summary:
# For each lab directory, we want to create a separate executable.
# e.g. create `build/lab2.exe` from `src/lab2/*.c` and `src/lab2/*/*.c` files.
#
# This makefile should (hopefully) work on both Windows and UNIX-like systems,
# without needing to set up a specific development environment (e.g. by using
# Docker). You will still need to have build tools installed and available on
# your system's shell (i.e. make, gcc, etc.)
#
# An alternative approach to take if this breaks would be to make a development
# environment with a POSIX `sh` shell using Docker, WSL, Cygwin, or similar.
# e.g. Dockerfile https://github.com/skeeto/w64devkit
# e.g. Open cygwin shell `cygwin`, then navigate to code directory on cygdrive
#      `cd C:/<path-on-windows>`.

# Additional makefile resources:
#
# RTFM
# https://www.gnu.org/software/make/manual/
#
# Readable makefile tutorial. No hints about Windows oddities.
# https://makefiletutorial.com/
#
# Here's a simple, documented example of a makefile for c projects on windows
# https://github.com/Goldan32/c-build-systems/tree/master/windows
# And some associated reddit discussion
# https://old.reddit.com/r/C_Programming/comments/uv10l2/makefile_example_project_for_windows_with_source/
#
# VEX Robotics have a cross-platform makefile build system for VEX V5 robots. 
# It's split across three files: makefile, mkenv.mk and mkrules.mk
# https://github.com/Wingus-Dingus-Robotics/WDR-Spin-Up/blob/main/makefile
#
# Portable makefiles... across POSIX systems. Did you know that GNU make isn't 
# the only implementation of make?
# https://nullprogram.com/blog/2017/08/20/

# Ensure that phony targets are not interpreted as file or directory names
.PHONY: all, clean, test, check, lab1, lab2, filter, project1, project2, project24, project25, project3

# DEBUG=1 -> Debug mode (no optimization, produce debug info)
# DEBUG=0 -> Release mode (full optimization, strip debug info)
DEBUG=1

# Choose shell for OS, and define shell command to remove directories
ifeq ($(OS),Windows_NT)
SHELL=cmd
RMDIR=rd /s /q
else
SHELL=sh
RMDIR=rm -rf
endif

# Choose C and C++ Compiler, linker, etc...
# e.g. gcc, clang, zig cc
#
# g++ is just gcc with some additional flags
# https://stackoverflow.com/questions/172587/what-is-the-difference-between-g-and-gcc
CC=gcc
CPP=g++
LINK=g++

# For a summary of gcc compiler options, RTFM.
# https://gcc.gnu.org/onlinedocs/gcc/Option-Summary.html

# C and C++ compiler options (warnings, iso standard)
CFLAGS += -Wall -Wextra -Wpedantic \
          -Wformat=2 -Wno-unused-parameter -Wshadow \
          -Wwrite-strings -Wstrict-prototypes -Wold-style-definition \
          -Wredundant-decls -Wnested-externs -Wmissing-include-dirs
CFLAGS += -std=c11
CPPFLAGS += -Wall -Wextra -Wpedantic \
			-Wformat=2 -Wno-unused-parameter -Wshadow \
			-Wwrite-strings -Wredundant-decls -Wmissing-include-dirs \
			-Wsign-conversion
CPPFLAGS += -std=c++14

# C and C++ compiler options (optimization level, debug)
ifneq ($(DEBUG), 0)
CFLAGS += -O0 -g
CPPFLAGS += -O0 -g
else
CFLAGS += -O2
CPPFLAGS += O2
endif

# Directories
BUILDDIR=build
SRCDIR=src

#################
# Build targets #
#################

# Target executable for each lab
# Executables should have `.exe` file extension on Windows.
ifeq ($(OS),Windows_NT)
LAB1:=$(BUILDDIR)/lab1/lab1.exe
LAB2:=$(BUILDDIR)/lab2/lab2.exe
FILTEREXAMPLE:=$(BUILDDIR)/filtering_example/filter.exe
PROJECT1:=$(BUILDDIR)/project1/project1.exe
PROJECT2:=$(BUILDDIR)/project2/project2.exe
PROJECT24:=$(BUILDDIR)/project24/project24.exe
PROJECT25:=$(BUILDDIR)/project25/project25.exe
PROJECT3:=$(BUILDDIR)/project3/project3.exe
else
LAB1:=$(BUILDDIR)/lab1/lab1
LAB2:=$(BUILDDIR)/lab2/lab2
FILTEREXAMPLE:=$(BUILDDIR)/filtering_example/filter
PROJECT1:=$(BUILDDIR)/project1/project1
PROJECT2:=$(BUILDDIR)/project2/project2
PROJECT24:=$(BUILDDIR)/project24/project24
PROJECT25:=$(BUILDDIR)/project25/project25
PROJECT3:=$(BUILDDIR)/project3/project3
endif

# C and C++ source and include files
# The "shaun_bmp" library should be built for all projects
SHAUNBMPSRC += $(wildcard $(SRCDIR)/shaun_bmp/*.c)
SHAUNBMPSRC += $(wildcard $(SRCDIR)/shaun_bmp/*/*.c)
SHAUNBMPSRC += $(wildcard $(SRCDIR)/shaun_bmp/*.cpp)
SHAUNBMPSRC += $(wildcard $(SRCDIR)/shaun_bmp/*/*.cpp)
SHAUNBMPINC += $(wildcard $(SRCDIR)/shaun_bmp/*.h)
SHAUNBMPINC += $(wildcard $(SRCDIR)/shaun_bmp/*/*.h)
SHAUNBMPINC += $(wildcard $(SRCDIR)/shaun_bmp/*.hpp)
SHAUNBMPINC += $(wildcard $(SRCDIR)/shaun_bmp/*/*.hpp)

LAB1SRC += $(wildcard $(SRCDIR)/lab1/*.c)
LAB1SRC += $(wildcard $(SRCDIR)/lab1/*/*.c)
LAB1SRC += $(wildcard $(SRCDIR)/lab1/*.cpp)
LAB1SRC += $(wildcard $(SRCDIR)/lab1/*/*.cpp)
LAB1INC += $(wildcard $(SRCDIR)/lab1/*.h)
LAB1INC += $(wildcard $(SRCDIR)/lab1/*/*.h)
LAB1INC += $(wildcard $(SRCDIR)/lab1/*.hpp)
LAB1INC += $(wildcard $(SRCDIR)/lab1/*/*.hpp)

LAB2SRC += $(wildcard $(SRCDIR)/lab2/*.c)
LAB2SRC += $(wildcard $(SRCDIR)/lab2/*/*.c)
LAB2INC += $(wildcard $(SRCDIR)/lab2/*.h)
LAB2INC += $(wildcard $(SRCDIR)/lab2/*/*.h)
LAB2SRC += $(wildcard $(SRCDIR)/lab2/*.cpp)
LAB2SRC += $(wildcard $(SRCDIR)/lab2/*/*.cpp)
LAB2INC += $(wildcard $(SRCDIR)/lab2/*.hpp)
LAB2INC += $(wildcard $(SRCDIR)/lab2/*/*.hpp)

FILTEREXAMPLESRC += $(wildcard $(SRCDIR)/filtering_example/*.cpp)
FILTEREXAMPLESRC += $(wildcard $(SRCDIR)/filtering_example/*/*.cpp)
FILTEREXAMPLEINC += $(wildcard $(SRCDIR)/filtering_example/*.h)
FILTEREXAMPLEINC += $(wildcard $(SRCDIR)/filtering_example/*/*.h)

PROJECT1SRC += $(wildcard $(SRCDIR)/project1/*.c)
PROJECT1SRC += $(wildcard $(SRCDIR)/project1/*/*.c)
PROJECT1SRC += $(wildcard $(SRCDIR)/project1/*.cpp)
PROJECT1SRC += $(wildcard $(SRCDIR)/project1/*/*.cpp)
PROJECT1INC += $(wildcard $(SRCDIR)/project1/*.h)
PROJECT1INC += $(wildcard $(SRCDIR)/project1/*/*.h)
PROJECT1INC += $(wildcard $(SRCDIR)/project1/*.hpp)
PROJECT1INC += $(wildcard $(SRCDIR)/project1/*/*.hpp)

PROJECT2SRC += $(wildcard $(SRCDIR)/project2/*.c)
PROJECT2SRC += $(wildcard $(SRCDIR)/project2/*/*.c)
PROJECT2SRC += $(wildcard $(SRCDIR)/project2/*.cpp)
PROJECT2SRC += $(wildcard $(SRCDIR)/project2/*/*.cpp)
PROJECT2INC += $(wildcard $(SRCDIR)/project2/*.h)
PROJECT2INC += $(wildcard $(SRCDIR)/project2/*/*.h)
PROJECT2INC += $(wildcard $(SRCDIR)/project2/*.hpp)
PROJECT2INC += $(wildcard $(SRCDIR)/project2/*/*.hpp)

PROJECT24SRC += $(wildcard $(SRCDIR)/project24/*.c)
PROJECT24SRC += $(wildcard $(SRCDIR)/project24/*/*.c)
PROJECT24SRC += $(wildcard $(SRCDIR)/project24/*.cpp)
PROJECT24SRC += $(wildcard $(SRCDIR)/project24/*/*.cpp)
PROJECT24INC += $(wildcard $(SRCDIR)/project24/*.h)
PROJECT24INC += $(wildcard $(SRCDIR)/project24/*/*.h)
PROJECT24INC += $(wildcard $(SRCDIR)/project24/*.hpp)
PROJECT24INC += $(wildcard $(SRCDIR)/project24/*/*.hpp)

PROJECT25SRC += $(wildcard $(SRCDIR)/project25/*.c)
PROJECT25SRC += $(wildcard $(SRCDIR)/project25/*/*.c)
PROJECT25SRC += $(wildcard $(SRCDIR)/project25/*.cpp)
PROJECT25SRC += $(wildcard $(SRCDIR)/project25/*/*.cpp)
PROJECT25INC += $(wildcard $(SRCDIR)/project25/*.h)
PROJECT25INC += $(wildcard $(SRCDIR)/project25/*/*.h)
PROJECT25INC += $(wildcard $(SRCDIR)/project25/*.hpp)
PROJECT25INC += $(wildcard $(SRCDIR)/project25/*/*.hpp)

PROJECT3SRC += $(wildcard $(SRCDIR)/project3/*.c)
PROJECT3SRC += $(wildcard $(SRCDIR)/project3/*/*.c)
PROJECT3SRC += $(wildcard $(SRCDIR)/project3/*.cpp)
PROJECT3SRC += $(wildcard $(SRCDIR)/project3/*/*.cpp)
PROJECT3INC += $(wildcard $(SRCDIR)/project3/*.h)
PROJECT3INC += $(wildcard $(SRCDIR)/project3/*/*.h)
PROJECT3INC += $(wildcard $(SRCDIR)/project3/*.hpp)
PROJECT3INC += $(wildcard $(SRCDIR)/project3/*/*.hpp)

# Object files to be built for each source file (both .c and .cpp)
# - Extract just the filenames without extensions from source filepaths
# - Add the desired build directory as filepath prefix (addprefix)
# - Add object file extension .o as file suffix (addsuffix)
# LAB1OBJ += $(addprefix $(BUILDDIR)/lab1/, $(notdir $(LAB1SRC:.c=.o)))
LAB1OBJ += 	$(addsuffix .o, \
			$(addprefix $(BUILDDIR)/lab1/, $(notdir $(basename $(LAB1SRC)))))
LAB2OBJ += 	$(addsuffix .o, \
			$(addprefix $(BUILDDIR)/lab2/, $(notdir $(basename $(LAB2SRC)))))

SHAUNBMPOBJ += 	$(addsuffix .o, \
				$(addprefix $(BUILDDIR)/shaun_bmp/, \
				$(notdir $(basename $(SHAUNBMPSRC)))))

FILTEREXAMPLEOBJ := $(addsuffix .o, \
					$(addprefix $(BUILDDIR)/filtering_example/, \
					$(notdir $(basename $(FILTEREXAMPLESRC)))))

PROJECT1OBJ += 	$(addsuffix .o, \
				$(addprefix $(BUILDDIR)/project1/, \
				$(notdir $(basename $(PROJECT1SRC)))))

PROJECT2OBJ += 	$(addsuffix .o, \
				$(addprefix $(BUILDDIR)/project2/, \
				$(notdir $(basename $(PROJECT2SRC)))))

PROJECT24OBJ += $(addsuffix .o, \
				$(addprefix $(BUILDDIR)/project24/, \
				$(notdir $(basename $(PROJECT24SRC)))))

PROJECT25OBJ += $(addsuffix .o, \
				$(addprefix $(BUILDDIR)/project25/, \
				$(notdir $(basename $(PROJECT25SRC)))))

PROJECT3OBJ += 	$(addsuffix .o, \
				$(addprefix $(BUILDDIR)/project3/, \
				$(notdir $(basename $(PROJECT3SRC)))))

#####################################################################################TODO finish this....

# The prerequisites to phony target "all" are the final executables to be built.
# Because "all" is the first target in the makefile, running the command `make`
# in a terminal will be equivalent to `make all`.
all: lab1 lab2 filter project1 project2 project24 project25 project3

# Build directory must exist as a prerequisite for every other target.
# However, we don't want to force every target to update whenever the contents
# of the build directory changes.
# To make `BUILDDIR` an order-only prerequisite, use the following rule format:
# ```<target>: <normal-prequisites> | $(BUILDDIR)```
$(BUILDDIR):
	mkdir $(BUILDDIR)
ifeq ($(OS),Windows_NT)
	mkdir $(BUILDDIR)\lab1
	mkdir $(BUILDDIR)\lab2
	mkdir $(BUILDDIR)\shaun_bmp
	mkdir $(BUILDDIR)\filtering_example
	mkdir $(BUILDDIR)\project1
	mkdir $(BUILDDIR)\project2
	mkdir $(BUILDDIR)\project24
	mkdir $(BUILDDIR)\project25
	mkdir $(BUILDDIR)\project3
else
	mkdir $(BUILDDIR)/lab1
	mkdir $(BUILDDIR)/lab2
	mkdir $(BUILDDIR)/shaun_bmp
	mkdir $(BUILDDIR)/filtering_example
	mkdir $(BUILDDIR)/project1
	mkdir $(BUILDDIR)/project2
	mkdir $(BUILDDIR)/project24
	mkdir $(BUILDDIR)/project25
	mkdir $(BUILDDIR)/project3
endif

# The build process for each target executable should be roughly the same:
# 1. Compile each source file into an object file separately, without linking.
#    This is achieved using the `-c` flag.
# 2. Link all of the object files together into the final executable.

# Build rules for shaun_bmp.o (library dependency for other projects)
shaun: $(SHAUNBMPOBJ)

$(BUILDDIR)/shaun_bmp/%.o: $(SRCDIR)/shaun_bmp/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/shaun_bmp/%.o: $(SRCDIR)/shaun_bmp/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/shaun_bmp/%.o: $(SRCDIR)/shaun_bmp/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/shaun_bmp/%.o: $(SRCDIR)/shaun_bmp/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for lab1.exe
lab1: $(LAB1)

$(LAB1): $(LAB1OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(LAB1OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/lab1/%.o: $(SRCDIR)/lab1/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/lab1/%.o: $(SRCDIR)/lab1/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/lab1/%.o: $(SRCDIR)/lab1/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/lab1/%.o: $(SRCDIR)/lab1/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for lab2.exe
lab2: $(LAB2)

$(LAB2): $(LAB2OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(LAB2OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/lab2/%.o: $(SRCDIR)/lab2/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/lab2/%.o: $(SRCDIR)/lab2/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/lab2/%.o: $(SRCDIR)/lab2/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/lab2/%.o: $(SRCDIR)/lab2/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for filter.exe (example code for lab2)
filter: $(FILTEREXAMPLE)

$(FILTEREXAMPLE): $(FILTEREXAMPLEOBJ)
	$(LINK) -o $@ $(FILTEREXAMPLEOBJ)

$(BUILDDIR)/filtering_example/%.o: $(SRCDIR)/filtering_example/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/filtering_example/%.o: $(SRCDIR)/filtering_example/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for project1.exe
project1: $(PROJECT1)

$(PROJECT1): $(PROJECT1OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(PROJECT1OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/project1/%.o: $(SRCDIR)/project1/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project1/%.o: $(SRCDIR)/project1/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project1/%.o: $(SRCDIR)/project1/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/project1/%.o: $(SRCDIR)/project1/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for project2.exe
project2: $(PROJECT2)

$(PROJECT2): $(PROJECT2OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(PROJECT2OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/project2/%.o: $(SRCDIR)/project2/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project2/%.o: $(SRCDIR)/project2/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project2/%.o: $(SRCDIR)/project2/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/project2/%.o: $(SRCDIR)/project2/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for project24.exe
project24: $(PROJECT24)

$(PROJECT24): $(PROJECT24OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(PROJECT24OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/project24/%.o: $(SRCDIR)/project24/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project24/%.o: $(SRCDIR)/project24/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project24/%.o: $(SRCDIR)/project24/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/project24/%.o: $(SRCDIR)/project24/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for project25.exe
project25: $(PROJECT25)

$(PROJECT25): $(PROJECT25OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(PROJECT25OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/project25/%.o: $(SRCDIR)/project25/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project25/%.o: $(SRCDIR)/project25/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project25/%.o: $(SRCDIR)/project25/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/project25/%.o: $(SRCDIR)/project25/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Build rules for project3.exe
project3: $(PROJECT3)

$(PROJECT3): $(PROJECT3OBJ) $(SHAUNBMPOBJ)
	$(LINK) -o $@ $(PROJECT3OBJ) $(SHAUNBMPOBJ)

$(BUILDDIR)/project3/%.o: $(SRCDIR)/project3/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project3/%.o: $(SRCDIR)/project3/*/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/project3/%.o: $(SRCDIR)/project3/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(BUILDDIR)/project3/%.o: $(SRCDIR)/project3/*/%.cpp | $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -c -o $@ $<

#####################
# Not build targets #
#####################

# The phony target "clean" will remove the build directory if it exists.
clean:
ifneq ("$(wildcard $(BUILDDIR))", "")
	$(RMDIR) $(BUILDDIR)
endif

# Detect bugs in source code using a static code analysis tool.
check:
	cppcheck src --check-level=exhaustive

# List dynamic dependencies required by compiled executables.
exe_list=$(LAB1) $(LAB2)
ldd:
ifeq ($(OS),Windows_NT)
	@$(foreach exe, $(exe_list), \
	echo $(exe) &\
	objdump -p $(exe) | FINDSTR /C:"DLL Name" &)
else
	@$(foreach exe, $(exe_list), \
	echo $(exe); \
	objdump -p $(exe) | grep NEEDED;)
endif
