#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shaun_bmp/bmp_io.h"

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: lab1 [options] <file>\n\
\n\
Takes a bitmap image file input <file> and writes a modified output file.\n\
\n\
If the -n option is specified, the a number is appended to the end of <file> \
for each sequential input.\n\
e.g. file.bmp -> file1.bmp, file2.bmp, etc.\n\
\n\
Option      Description\n\
-o <file>   Specify output file (default: output.bmp)\n\
-n <int>    Number of sequential input files";

int main(int argc, char *argv[]) {
    // Print help if no arguments provided
    if (argc < 1) {
        puts(CLI_HELP);
        return EXIT_FAILURE;
    }

    // Default command options
    char * inputFile = NULL;    // <inputFile>
    char * outputFile = NULL;   // -o <outputFile>
    int numberOfInputs = 0;     // -n <number>

    // Parse command-line arguments
    enum { O, N, NOT_OPTION } currentArg, prevArg = NOT_OPTION;
    for (int i_arg = 1; i_arg < argc; i_arg++) {
        // Check if current argument is an option flag. If so, then the next 
        // argument will contain the value of that option.
        if (strcmp(argv[i_arg], "-o") == 0) {
            currentArg = O;
        }
        else if (strcmp(argv[i_arg], "-n") == 0) {
            currentArg = N;
        }
        else {
            currentArg = NOT_OPTION;
        }

        // Set option value (if preceding argument was an option flag).
        if (prevArg == O) {
            outputFile = argv[i_arg];
        }
        else if (prevArg == N) {
            numberOfInputs = atoi(argv[i_arg]);
        }

        // Check if current argument is the input file. It could be the first 
        // or last argument. Last takes precedence, since that's the format on 
        // the CLI help message.
        //
        // Case 1:
        // If the first argument is not a flag, then it must be the input file.
        //
        // Case 2: 
        // If the final argument is not an option value or a flag, then it must 
        // be the input file.
        if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
            inputFile = argv[i_arg];
        }
        if ((i_arg == (argc - 1)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
            inputFile = argv[i_arg];
        }

        // Store state of current argument
        prevArg = currentArg;
    }

    // Print help if no input file was provided
    if (inputFile == NULL) {
        puts(CLI_HELP);
        return EXIT_FAILURE;
    }
    
    // Hey look, it works!
    printf("i = %s\n", inputFile);
    printf("n = %d\n", numberOfInputs);
    printf("o = %s\n", outputFile);

    // TODO: Run the actual application
    int fileErr = read_header(inputFile);
    if (fileErr == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
    else if (fileErr == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
    else if (fileErr == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
    else if (fileErr == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
    else if (fileErr == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");


    return EXIT_SUCCESS;
}
