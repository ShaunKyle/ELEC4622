#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "../shaun_bmp/bmp_io.h"

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: project3 [options] <source> <target> <B> <S>\n\
\n\
Estimates motion between the <source> and <target> frames.\n\
Assumes that entire blocks of pixels undergo the same motion.\n\
The output image is the target frame annotated with motion vectors.\n\
\n\
Parameters:\n\
<B> specifies the block size.\n\
<S> specifies the motion search range [-S,+S] in each direction.\n\
\n\
Option      Description\n\
-o <file>   Specify output file (default: output.bmp)\n\
-m          Change motion vector search criterion to minimize MSE instead of \
SAD\n\
-w          Change interpolation method to windowed sinc instead of bilinear\n\
";

int main (int argc, char *argv[]) {
    /////////////////////
    // CLI boilerplate //
    /////////////////////

    // Print help if no arguments provided
    if (argc < 1) {
        puts(CLI_HELP);
        return EXIT_FAILURE;
    }

    // Default command options
    char * sourceFile = NULL;   // <source frame>
    char * targetFile = NULL;   // <target frame>
    int B = 0;                  // <B>
    int S = 0;                  // <S>
    char * outputFile = NULL;                   // -o <outputFile>
    bool minimizeMSEFlag = false;               // -m
    bool interpolateWindowedSincFlag = false;   // -w

    // Parse command-line arguments
    enum { O, M, W, NOT_OPTION } currentArg, prevArg = NOT_OPTION;
    for (int i_arg = 1; i_arg < argc; i_arg++) {
        // Check if current argument is an option flag. If so, then the next 
        // argument will contain the value of that option.
        if (strcmp(argv[i_arg], "-o") == 0) {
            currentArg = O;
        }
        else if (strcmp(argv[i_arg], "-m") == 0) {
            currentArg = M;
            minimizeMSEFlag = true;
        }
        else if (strcmp(argv[i_arg], "-w") == 0) {
            currentArg = W;
            interpolateWindowedSincFlag = true;
        }
        else {
            currentArg = NOT_OPTION;
        }

        // Set option value (if preceding argument was an option flag).
        if (prevArg == O) {
            outputFile = argv[i_arg];
        }
        // else if (prevArg == N) {
        //     numberOfInputs = atoi(argv[i_arg]);
        // }

        // Check if current argument is the <source> file. Could be the first 
        // or last argument.
        //
        // Case 1:
        // If the first argument is not a flag, then it must be <source>.
        // The second argument will be <target>.
        // The third argument will be <B>.
        // The fourth argument will be <S>.
        //
        // Case 2: 
        // If the 4th-to-last argument is not an option value or a flag, then 
        // it must be <source>.
        if (sourceFile == NULL) {
            if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
                // puts("case 1");
                sourceFile = argv[i_arg];
                targetFile = argv[i_arg+1];
                B = atoi(argv[i_arg+2]);
                S = atoi(argv[i_arg+3]);
            }
            else if ((i_arg == (argc - 4)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
                // puts("case 2");
                sourceFile = argv[i_arg];
                targetFile = argv[i_arg+1];
                B = atoi(argv[i_arg+2]);
                S = atoi(argv[i_arg+3]);
            }
        }

        // Store state of current argument
        prevArg = currentArg;
    }

    // Print help if no input files are provided
    if ((sourceFile == NULL) || (targetFile == NULL)) {
        puts(CLI_HELP);
        return EXIT_FAILURE;
    }

    // Give output file a default name if none provided
    if (outputFile == NULL) {
        outputFile = (char *) "out.bmp";
    }
    
    // Hey look, it works!
    puts("Input and output files:");
    printf("  <source> = %s\n", sourceFile);
    printf("  <target> = %s\n", targetFile);
    printf("  -o = %s\n\n", outputFile);

    puts("Parameters and options:");
    printf("  <B> = %d\n", B);
    printf("  <S> = %d\n", S);
    printf("  -m = %s\n\n", minimizeMSEFlag ? "MSE" : "SAD");
    printf("  -w = %s\n\n", interpolateWindowedSincFlag ? 
        "Windowed sinc interpolator" : "Bilinear interpolator");

    // Load input bitmap files
    bmp source_bmp, target_bmp;
    int fileErr = load_bmp(&source_bmp, sourceFile);
    if (fileErr != 0) {
        print_bmp_file_error(fileErr);
        return EXIT_FAILURE;
    }
    fileErr = load_bmp(&target_bmp, targetFile);
    if (fileErr != 0) {
        print_bmp_file_error(fileErr);
        return EXIT_FAILURE;
    }
    
    // printf("Info about %s\n", inputFile);
    // printf("Width x Height: %dx%d px\n", input_bmp.cols, input_bmp.rows);
    // printf("Components: %d\n", input_bmp.num_components);
    // printf("Line bytes: %d\n", input_bmp.line_bytes);
    // printf("Alignment padding bytes: %d\n\n", input_bmp.alignment_bytes);

    //////////////////////////////////////
    // Task 1: Calculate motion vectors //
    //////////////////////////////////////

    // SAD = sum of absolute differences
    // MSE = mean square error



    return EXIT_SUCCESS;
}
