#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "../shaun_bmp/bmp_io.h"

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: project2 [options] <file> <D> <H>\n\
\n\
Creates a Laplacian domain representation of an image internally, then \n\
reconstructs the original image from the Laplacian domain representation. \n\
\n\
Create Laplacian pyramid:\n\
- Writes out <D>+1 sub-images (x0, x1, ..., xD) as a single BMP.\n\
- Region of support for windowed sincs (interpolation kernel) is [-H, H].\n\
\n\
Reconstruct image from Laplacian pyramid:\n\
- \n\
\n\
Option      Description\n\
-o <file>   Specify output file (default: output.bmp)\n\
-g          Gaussian pyramid instead of Laplacian.\n\
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
    char * inputFile = NULL;    // <inputFile>
    int D = 0;                  // <D>
    int H = 0;                  // <H>
    char * outputFile = NULL;           // -o <outputFile>
    bool gaussianPyramidFlag = false;   // -g

    // Parse command-line arguments
    enum { O, G, NOT_OPTION } currentArg, prevArg = NOT_OPTION;
    for (int i_arg = 1; i_arg < argc; i_arg++) {
        // Check if current argument is an option flag. If so, then the next 
        // argument will contain the value of that option.
        if (strcmp(argv[i_arg], "-o") == 0) {
            currentArg = O;
        }
        else if (strcmp(argv[i_arg], "-g") == 0) {
            currentArg = G;
            gaussianPyramidFlag = true;
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

        // Check if current argument is the input file. It could be the first 
        // or last argument. Last takes precedence, since that's the format on 
        // the CLI help message.
        //
        // Case 1:
        // If the first argument is not a flag, then it must be the input file.
        // The second argument will be <D>.
        // The third argument will be <H>.
        //
        // Case 2: 
        // If the 3rd-to-last argument is not an option value or a flag, then 
        // it must be the input file (final args will be <D> and <H>).
        if (inputFile == NULL) {
            if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
                // puts("case 1");
                inputFile = argv[i_arg];
                D = atoi(argv[i_arg+1]);
                H = atoi(argv[i_arg+2]);
            }
            else if ((i_arg == (argc - 3)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
                // puts("case 2");
                inputFile = argv[i_arg];
                D = atoi(argv[i_arg+1]);
                H = atoi(argv[i_arg+2]);
            }
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
    printf("o = %s\n\n", outputFile);

    printf("D = %d\n", D);
    printf("H = %d\n", H);
    printf("G = %s\n\n", gaussianPyramidFlag ? "Gaussian" : "Laplacian");

    // Load input bitmap file
    bmp input_bmp;
    int fileErr = load_bmp(&input_bmp, inputFile);
    if (fileErr != 0) {
        print_bmp_file_error(fileErr);
        return EXIT_FAILURE;
    }
    
    printf("Info about %s\n", inputFile);
    printf("Width x Height: %dx%d px\n", input_bmp.cols, input_bmp.rows);
    printf("Components: %d\n", input_bmp.num_components);
    printf("Line bytes: %d\n", input_bmp.line_bytes);
    printf("Alignment padding bytes: %d\n\n", input_bmp.alignment_bytes);


    ///////////////////////////////////////
    // Task 1: Create a Gaussian Pyramid //
    ///////////////////////////////////////

    if (gaussianPyramidFlag) {
        puts("TODO> Task 1");
    }


    ////////////////////////////////////////
    // Task 2: Create a Laplacian Pyramid //
    ////////////////////////////////////////

    if (!gaussianPyramidFlag) {
        puts("TODO> Task 2");
    }


    ////////////////////////////
    // Output for task 1 or 2 //
    ////////////////////////////
    

    return EXIT_SUCCESS;
}
