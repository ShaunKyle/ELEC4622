#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h> 

#include "../shaun_bmp/bmp_io.h"
#include "../shaun_bmp/image.h"

float hanning_window(int n, int extent);
void design_windowed_sinc_filter(pixel_t *h_sinc, int extent, float stretch);

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
        puts("Task 1> Gaussian Pyramid");        

        // Design windowed stretched sinc
        // Reminder: Use apply_separable_filters_2n() to apply filter
        pixel_t h_sinc[2*H+1];
        design_windowed_sinc_filter(h_sinc, H, 2.0);

        // Allocate memory for Gaussian pyramid, based on height of input image
        image imgGaussianPyramid, image0;
        read_image_from_bmp(&image0, &input_bmp, H);
        perform_boundary_extension(&image0);
        int total_height = 0;
        for (int level = 0; level <= D; level++) {
            total_height += (image0.rows / pow(2, level));
        }
        printf("Total height of pyramid: %d\n", total_height);
        imgGaussianPyramid.rows = total_height;
        imgGaussianPyramid.cols = image0.cols;
        imgGaussianPyramid.border = 0;
        imgGaussianPyramid.stride = image0.cols;
        imgGaussianPyramid.num_components = image0.num_components;
        imgGaussianPyramid.handle = malloc(
            total_height * imgGaussianPyramid.cols * 
            imgGaussianPyramid.num_components * sizeof(pixel_t)
        );
        imgGaussianPyramid.buf = imgGaussianPyramid.handle;

        // Initialize pyramid pixels to 0
        for (int row = 0; row < imgGaussianPyramid.rows; row++) {
            for (int col = 0; col < imgGaussianPyramid.cols; col++) {
                const int index = 
                    (row * imgGaussianPyramid.stride + col) 
                    * imgGaussianPyramid.num_components;
                for (int plane = 0; plane < imgGaussianPyramid.num_components; plane++) {
                    imgGaussianPyramid.buf[index+plane] = 0.0;
                }
            }
        }

        // Create sub-images in pyramid
        int accessed_rows = 0;
        image imagePrev, imageCurrent;
        imageCurrent = image0;
        for (int level = 0; level <= D; level++) {
            printf("Level %d\n", level);
            printf("Starting row: %d\n", accessed_rows);

            // Write sub-image to pyramid
            for (int row = 0; row < imageCurrent.rows; row++) {
                for (int col = 0; col < imageCurrent.cols; col++) {
                    // Pixel to read from imageCurrent
                    const int row_flipped = imageCurrent.rows - row - 1;
                    const int indexCurrent = 
                        (row_flipped * imageCurrent.stride + col) 
                        * imageCurrent.num_components;
                    
                    // Calculate pyramid coordinates and array index
                    // Flip n2 because of BMP coordinate system
                    const int n1 = col;
                    const int n2 = total_height - accessed_rows - 1;
                    const int indexPyramid = 
                        (n2 * imgGaussianPyramid.stride + n1) 
                        * imgGaussianPyramid.num_components;
                    
                    // Write pixel from imageCurrent to pyramid
                    for (int plane = 0; plane < imageCurrent.num_components; plane++) {
                        imgGaussianPyramid.buf[indexPyramid+plane] = 
                            imageCurrent.buf[indexCurrent+plane];
                    }
                }

                // Ready to write next row
                // printf("%d\n", accessed_rows);
                accessed_rows++;
            }

            // Apply decimation filter to image (i.e. every second pixel)
            copy_image(&imageCurrent, &imagePrev, H);
            perform_boundary_extension(&imagePrev);
            apply_separable_filters_2n(&imagePrev, &imageCurrent, 
                h_sinc, h_sinc, H, H
            );
        }

        // Output
        export_image_as_bmp(&imgGaussianPyramid, outputFile);
        
    }


    ////////////////////////////////////////
    // Task 2: Create a Laplacian Pyramid //
    ////////////////////////////////////////

    if (!gaussianPyramidFlag) {
        puts("TODO> Task 2");

        // Design windowed sinc interpolator
        // Confine bandwidth to (-pi/2, pi/2)^2
        pixel_t g_interp[3];

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


void design_windowed_sinc_filter(pixel_t *h_sinc, int extent, float stretch) {
    // Design windowed sinc
    const int SINC_DIM = (2*extent+1);
    for (int tap = 0; tap < SINC_DIM; tap++) {
        const int n = tap - extent;
        if (n == 0) {
            h_sinc[tap] = 1;
        } else {
            h_sinc[tap] = sin(M_PI * n / stretch) / (M_PI * n / stretch);
        }
        h_sinc[tap] *= hanning_window(n, extent);
        printf("h_sinc[%d] = %f\n", n, h_sinc[tap]);
    }

    // Check DC gain of windowed sinc
    double dc_gain = 0;
    for (int row = 0; row < SINC_DIM; row++) {
        for (int col = 0; col < SINC_DIM; col++) {
            dc_gain += (h_sinc[row] * h_sinc[col]);
        }
    }
    printf("DC gain: %f\n", dc_gain);

    // Normalize dc_gain to 1
    // Note that this design is for separable filters. Sqrt.
    for (int tap = 0; tap < SINC_DIM; tap++) {
        h_sinc[tap] /= sqrt(dc_gain);
    }

    // Check normalized DC gain of windowed sinc
    double dc_gain_normalized = 0;
    for (int row = 0; row < SINC_DIM; row++) {
        for (int col = 0; col < SINC_DIM; col++) {
            dc_gain_normalized += (h_sinc[row] * h_sinc[col]);
        }
    }
    printf("DC gain normalized: %f\n", dc_gain_normalized);
}
