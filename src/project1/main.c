#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>   // fabs

#include "../shaun_bmp/bmp_io.h"
#include "../shaun_bmp/image.h"

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: project1 [options] <file> <sigma>\n\
\n\
Sigma is a floating point value. Default filter is Gaussian.\n\
\n\
Option      Description\n\
-o <file>   Specify output file (default: output.bmp)\n\
-w          Switch to moving average filter (default: Gaussian)\n\
-n <int>    Number of sequential input files";

int main(int argc, char *argv[]) {
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
    float sigma = 0.0;          // <sigma>
    char * outputFile = NULL;       // -o <outputFile>
    bool movingAverageFlag = false; // -w
    int numberOfInputs = 0;         // -n <number>

    // Parse command-line arguments
    enum { O, N, W, NOT_OPTION } currentArg, prevArg = NOT_OPTION;
    for (int i_arg = 1; i_arg < argc; i_arg++) {
        // Check if current argument is an option flag. If so, then the next 
        // argument will contain the value of that option.
        if (strcmp(argv[i_arg], "-o") == 0) {
            currentArg = O;
        }
        else if (strcmp(argv[i_arg], "-n") == 0) {
            currentArg = N;
        }
        else if (strcmp(argv[i_arg], "-w") == 0) {
            currentArg = W;
            movingAverageFlag = true;
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
        // The second argument will be sigma.
        //
        // Case 2: 
        // If the 2nd to last argument is not an option value or a flag, then 
        // it must be the input file (and final argument will be sigma).
        if (inputFile == NULL) {
            if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
                puts("case 1");
                inputFile = argv[i_arg];
                sigma = atof(argv[i_arg+1]);
            }
            else if ((i_arg == (argc - 2)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
                puts("case 1");
                inputFile = argv[i_arg];
                sigma = atof(argv[i_arg+1]);
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

    printf("sigma = %f\n", sigma);
    printf("w = %s\n\n", movingAverageFlag ? "Moving Average" : "Gaussian");

    printf("n = %d\n\n", numberOfInputs);

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

    ///////////////////////////
    // More interesting part //
    ///////////////////////////

    if (movingAverageFlag) {
        // Moving Average filter with DC gain of 1

        // Determine appropriate region of support... using brute force.
        const float target_var = sigma * sigma;
        
        float abs_error = 10000000000000.0; // ridiculously implausible value
        float prev_abs_error = 0.0;

        int H = 0;  // Choose region of support [-H, +H]^2
        int DIM;
        int TAPS;
        while (true) {
            DIM = (2*H+1);
            TAPS = DIM * DIM;
            float current_var = 0.0;
            for (int n1 = -H; n1 < H+1; n1++) {
                for (int n2 = -H; n2 < H+1; n2++) {
                    current_var += (n1*n1 + n2*n2) * (1.0 / TAPS);
                }
            }
            prev_abs_error = abs_error;
            abs_error = fabs(target_var - current_var);
            printf("H=%d has Var error of %f\n", H, abs_error);

            if (abs_error > prev_abs_error) {
                H--;
                break;
            } else {
                H++;
            }
        }
        printf("Chose H=%d\n\n", H);

        DIM = (2*H+1);
        TAPS = DIM * DIM;
        

        // // Direct method
        // pixel_t psf_ma[TAPS];
        // float dc_gain = 0;
        // for (int row = 0; row < DIM; row++) {
        //     for (int col = 0; col < DIM; col++) {
        //         psf_ma[col + row*DIM] = 1.0 / TAPS;
        //         // printf("h[%d, %d] = %f\n", col, row, 1.0 / TAPS);
        //         dc_gain += psf_ma[col + row*DIM];
        //     }
        // }
        // printf("DC gain is %f\n", dc_gain);
        // printf("Region of support is [-%d, %d]^2\n", H, H);

        // // Use moving average filter
        // image imageIn, imageOut;
        // read_image_from_bmp(&imageIn, &input_bmp, H);
        // perform_boundary_extension(&imageIn);
        // apply_filter(&imageIn, &imageOut, psf_ma, H, H);
        // export_image_as_bmp(&imageOut, outputFile);

        // Bonus-1: Separable method
        pixel_t h1[DIM];
        pixel_t h2[DIM];
        float dc_gain_h1 = 0;
        float dc_gain = 0;
        for (int col = 0; col < DIM; col++) {
            h1[col] = 1.0 / TAPS;
            dc_gain_h1 += h1[col];
        }
        for (int row = 0; row < DIM; row++) {
            h2[row] = 1;
            dc_gain += dc_gain_h1 * h2[row];
        }

        printf("h1 region of support is [-%d, %d]^2\n", H, 1);
        printf("h2 region of support is [-%d, %d]^2\n", 1, H);
        printf("Total DC gain is %f\n", dc_gain);

        // Use cascaded separable moving average filters
        image imageIn, imageOut;
        read_image_from_bmp(&imageIn, &input_bmp, H);
        perform_boundary_extension(&imageIn);
        apply_separable_filters(&imageIn, &imageOut, h1, h2, H, H);
        export_image_as_bmp(&imageOut, outputFile);
        
    }
    else {
        // Design Gaussian filter
        // Some reading https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
        // Also notes Ch2 Section 4.1 Page 24
        const int H = 3*sigma;   // Choose region of support [-H, +H]^2
        const int DIM = (2*H+1);
        const int TAPS = DIM * DIM;
        pixel_t psf_gaussian[TAPS];

        // Direct calculation of Gaussian PSF values
        float dc_gain = 1.0;
        for (int row = 0; row < DIM; row++) {
            for (int col = 0; col < DIM; col++) {
                const int n1 = col - H;
                const int n2 = row - H;
                const float scale = 1.0 / (2 * M_PI * sigma * sigma);
                psf_gaussian[col + row*DIM] = scale * exp(
                    -(n1*n1 + n2*n2) / (2*sigma*sigma)
                );
                // printf("h[%d, %d] = %f\n", col, row, 1.0 / TAPS);
                dc_gain += psf_gaussian[col + row*DIM];
            }
        }
        printf("DC gain is %f\n", dc_gain);
        printf("Region of support is [-%d, %d]^2\n", H, H);

        // Use Gaussian filter
        image imageIn, imageOut;
        read_image_from_bmp(&imageIn, &input_bmp, H);
        perform_boundary_extension(&imageIn);
        apply_filter(&imageIn, &imageOut, psf_gaussian, H, H);
        export_image_as_bmp(&imageOut, outputFile);
    }

    // TODO: Handle -n flag.

    return EXIT_SUCCESS;
}
