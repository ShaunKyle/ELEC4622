# Multimedia Signal Processing

Practical implementation of signal processing algorithms in mostly C++, for the 
course ELEC4622 "Multimedia Signal Processing" at the University of New South 
Wales.

| Application | Summary                                                       |
|-------------|---------------------------------------------------------------|
| lab1        | Basic read and write operations for `.bmp` images, in C.      |
| lab2        | Implementing 2D FIR filters. SSE2-family vector instructions. |

## Quickstart

See [SETUP.md](./SETUP.md) for developer setup instructions. You will need to 
have the required programs installed.

Compile applications for each lab. This should work on Windows/Linux/MacOS.

```sh
# Compile every application
make all

# Compile only the application for lab 1
make lab1

# Optional: Generate a JSON compilation database for clangd
compiledb make all
```

Run CLI application from command-line. CLI help message should probably be 
displayed.

```sh
# Run lab 1 (Windows)
.\build\lab1\lab1.exe

# Run lab 1 (UNIX-like)
./build/lab1/lab1
```
