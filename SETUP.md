# Developer setup

## Required software

Have the following programs available on your PATH:

- make (build tool)
- gcc (c compiler)
- gdb (debugger)
- cppcheck (static code analysis)

On Windows, you can use a package manager such as [Scoop](https://scoop.sh/) 
to handle the installation.

I like using Scoop because it keeps all installed packages isolated in a single 
directory `~\scoop`. Other package managers for Windows include winget and 
Chocolatey.

```sh
scoop install make gcc gdb cppcheck
```

## Additional software

A gdb frontend like [gdbgui](https://www.gdbgui.com/) could be nice.

Having a method for generating a Clang JSON Compilation Database can be useful 
for setting up a C/C++ language server such as 
[clangd](https://clangd.llvm.org/) for your IDE. Some build systems such as 
CMake support generating this database natively. For a makefile project, you 
can use a tool such as bear or 
[compiledb](https://github.com/nickdiego/compiledb) to convert a build log into 
a compilation database.

```sh
pipx install compiledb
```

## Setting up VSCode as a C/C++ development environment

All you need is a text editor and a terminal. And maybe a vertical ruler to 
help enforce an 80 character line limit that keeps source files readable.

JKJK, here's a list of 999+ extensions to download and the latest AI coding 
assistant you'll need a subscription for:

- [clangd](https://marketplace.visualstudio.com/items?itemName=llvm-vs-code-extensions.vscode-clangd)
- [Native Debug](https://open-vsx.org/extension/webfreak/debug)

## Setup for UNIX-like OS

You'll figure it out.
