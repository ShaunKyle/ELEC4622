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

    // Give output file a default name if none provided
    if (outputFile == NULL) {
        outputFile = (char *) "out.bmp";
    }
    
    // Hey look, it works!
    printf("i = %s\n", inputFile);
    printf("n = %d\n", numberOfInputs);
    printf("o = %s\n\n", outputFile);

    // Load input bitmap file
    bmp input_bmp;
    int fileErr = load_bmp(&input_bmp, inputFile);
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
    
    printf("Info about %s\n", inputFile);
    printf("Width x Height: %dx%d px\n", input_bmp.cols, input_bmp.rows);
    printf("Components: %d\n", input_bmp.num_components);
    printf("Line bytes: %d\n", input_bmp.line_bytes);
    printf("Alignment padding bytes: %d\n", input_bmp.alignment_bytes);

    // Create new bitmap file and write contents of inputFile into it, but make 
    // pixels half as bright.
    bmp output_bmp;
    create_bmp(&output_bmp, outputFile, 
        input_bmp.cols, input_bmp.rows, input_bmp.num_components);
    
    while(input_bmp.num_unaccessed_rows > 0) {
        uint8_t *current_line = malloc(input_bmp.line_bytes);
        int fileReadErr = read_bmp_line(&input_bmp, current_line);
        if (fileReadErr != 0) {
            // TODO: Handle the error...
            puts("AAAAAAAAAAAAAAAARRRRRRRRRGGGGGGHHHH!!!!!!");
            printf("%d\n", fileReadErr);
            return EXIT_FAILURE;
        }
        // TODO: Write line to output_bmp...
    }

    close_bmp(&input_bmp);
    close_bmp(&output_bmp);

    return EXIT_SUCCESS;
}
