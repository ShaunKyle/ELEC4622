#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../shaun_bmp/bmp_io.h"

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: lab1 [options] <file>\n\
\n\
Takes a bitmap image file input <file> and writes a modified output file.\n\
\n\
If the -n option is specified, a number is appended to the end of <file> \
for each sequential input and corresponding output file.\n\
For example, the command `lab1 -n 2 file.bmp` will:\n\
- Read input files {file1.bmp, file2.bmp} \n\
- Write output files {out1.bmp, out2.bmp} \n\
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
    if (fileErr != 0) {
        print_bmp_file_error(fileErr);
        return EXIT_FAILURE;
    }
    
    printf("Info about %s\n", inputFile);
    printf("Width x Height: %dx%d px\n", input_bmp.cols, input_bmp.rows);
    printf("Components: %d\n", input_bmp.num_components);
    printf("Line bytes: %d\n", input_bmp.line_bytes);
    printf("Alignment padding bytes: %d\n", input_bmp.alignment_bytes);

    // What's the goal for this program? idk.
    //
    // Create new bitmap file and write contents of inputFile into it, but make 
    // some pixels half as bright, and some zero.

    // Read every line of input image into imageStore.data
    image imageStore;
    int readErr = read_image_from_bmp(&imageStore, &input_bmp);
    if (readErr != 0) {
        print_bmp_file_error(readErr);
        return EXIT_FAILURE;
    }

    // Modify each component of pixel data separately.
    const int width = imageStore.cols;
    const int height = imageStore.rows;
    const int planes = imageStore.num_components;
    uint8_t *image_data = imageStore.data;
    for (int nRow = 0; nRow < height; nRow++) {
        for (int nCol = 0; nCol < width; nCol++) {
            int n = nRow*(width * planes) + nCol*planes;
            if (planes == 3) {
                image_data[n] = (uint8_t) (image_data[n+2] * 0.5);      // blue
                image_data[n+1] = 0;    // green
                image_data[n+2] = 0;    // red
            }
            else {
                image_data[n] = (uint8_t) (image_data[n] * 0.5); // grayscale
            }
        }
    }

    // Export processed pixel data to a new bitmap file
    int writeErr = export_image_as_bmp(&imageStore, outputFile);
    if (writeErr != 0) {
        print_bmp_file_error(writeErr);
        return EXIT_FAILURE;
    }

    // TODO: Process -n flag.

    return EXIT_SUCCESS;
}
