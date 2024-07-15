#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h> 

#include "../shaun_bmp/bmp_io.h"
#include "../shaun_bmp/image.h"

float hanning_window(int n, int extent);

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

        image image0, image1;
        read_image_from_bmp(&image0, &input_bmp, H);
        perform_boundary_extension(&image0);

        // Low pass filter and decimate image using stretched sinc


        // Design windowed stretched sinc
        #define SINC_DIM    (2*H+1)
        pixel_t h_sinc[SINC_DIM];
        for (int tap = 0; tap < SINC_DIM; tap++) {
            const int n = tap - H;
            if (n == 0) {
                h_sinc[tap] = 1;
            } else {
                h_sinc[tap] = sin(M_PI * n / 2.0) / (M_PI * n/2.0);
            }
            h_sinc[tap] *= hanning_window(n, H);
            printf("h_sinc[%d] = %f\n", n, h_sinc[tap]);
        }

        // Check DC gain of windowed sinc
        double dc_gain = 0;
        for (int row = 0; row < SINC_DIM; row++) {
            for (int col = 0; col < SINC_DIM; col++) {
                // const int n1 = row - H;
                // const int n2 = col - H;
                dc_gain += (h_sinc[row] * h_sinc[col]);
            }
        }
        printf("DC gain: %f\n", dc_gain);

        // Normalize dc_gain to 1
        for (int tap = 0; tap < SINC_DIM; tap++) {
            h_sinc[tap] /= dc_gain;
        }

        // Apply filter, but to every second pixel.
        apply_separable_filters_2n(&image0, &image1, h_sinc, h_sinc, H, H);
        export_image_as_bmp(&image1, outputFile);

        
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


float hanning_window(int n, int extent) {
    // If extent is 0, we don't need a window.
    if (extent == 0) {
        return 1.0;
    }

    // Hanning window
    if (abs(n) < extent) {
        // printf("%f\n", cos(M_PI * n / extent));
        return (1 + cos(M_PI * n / extent)) / 2.0;
    }
    else {
        // puts("wall");
        return 0.0;
    }
}
