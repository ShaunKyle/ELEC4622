#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>   // fabs, exp, M_PI

#include "../shaun_bmp/bmp_io.h"
#include "../shaun_bmp/image.h"

// Debug flags
#define EXPORT_INTERMEDIATE_STEPS
// #define SHOW_TAP_VALUES

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: project1 [options] <file> <sigma> <alpha>\n\
\n\
Performs low pass filtering and differentiation on image. \n\
Sigma is std dev between 1.0 to 100.0\n\
Alpha is a real-valued positive scaling factor.\n\
\n\
Default low pass filter is Gaussian.\n\
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
    float alpha = 0.0;          // <alpha>
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
        // The third argument will be alpha.
        //
        // Case 2: 
        // If the 3rd-to-last argument is not an option value or a flag, then 
        // it must be the input file (and final args will be sigma and alpha).
        if (inputFile == NULL) {
            if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
                puts("case 1");
                inputFile = argv[i_arg];
                sigma = atof(argv[i_arg+1]);
                alpha = atof(argv[i_arg+2]);
            }
            else if ((i_arg == (argc - 3)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
                puts("case 1");
                inputFile = argv[i_arg];
                sigma = atof(argv[i_arg+1]);
                alpha = atof(argv[i_arg+2]);
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
    printf("alpha = %f\n", alpha);
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

    /////////////////////////////
    // Task 1: Low pass filter //
    /////////////////////////////

    //            +----------+
    // imageIn -->| Low pass |--> imageLowPass
    //            +----------+

    image imageIn, imageLowPass;

    if (movingAverageFlag) {
        // Moving Average filter with DC gain of 1
        puts("Task 1: Low pass filter (moving average)");

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
        read_image_from_bmp(&imageIn, &input_bmp, H);
        perform_boundary_extension(&imageIn);
        apply_separable_filters(&imageIn, &imageLowPass, h1, h2, H, H);
        
    }
    else {
        puts("Task 1: Low pass filter (Gaussian)");

        // Design Gaussian filter
        // Some reading https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
        // Also notes Ch2 Section 4.1 Page 24
        // Note that sampling period for digital image processing is T = 1
        const int H = 3*sigma;   // Choose region of support [-H, +H]^2
        const int DIM = (2*H+1);
        // const int TAPS = DIM * DIM;

        // // Direct calculation of 2D Gaussian PSF values
        // pixel_t psf_gaussian[TAPS];
        // float dc_gain = 1.0;
        // for (int row = 0; row < DIM; row++) {
        //     for (int col = 0; col < DIM; col++) {
        //         const int n1 = col - H;
        //         const int n2 = row - H;
        //         const float scale = 1.0 / (2 * M_PI * sigma * sigma);
        //         psf_gaussian[col + row*DIM] = scale * exp(
        //             -(n1*n1 + n2*n2) / (2*sigma*sigma)
        //         );
        //         printf("h[%d, %d] = %f\n", n1, n2, psf_gaussian[col+row*DIM]);
        //         dc_gain += psf_gaussian[col + row*DIM];
        //     }
        // }
        // printf("DC gain is %f\n", dc_gain);
        // printf("Region of support is [-%d, %d]^2\n", H, H);

        // // Use Gaussian filter
        // image imageIn, imageOut;
        // read_image_from_bmp(&imageIn, &input_bmp, H);
        // perform_boundary_extension(&imageIn);
        // apply_filter(&imageIn, &imageOut, psf_gaussian, H, H);

        // Bonus-1: Separable method
        pixel_t h_gaussian_1D[DIM];
        float dc_gain_1D = 0.0;
        for (int i = 0; i < DIM; i++) {
            const int n = i - H;
            const float scale = 1.0 / (2 * M_PI * sigma * sigma);
            h_gaussian_1D[i] = sqrt(scale) * exp(-n*n / (2*sigma*sigma));
            // Note: we use sqrt because filter is applied twice.

            #ifdef SHOW_TAP_VALUES
            printf("h[%d] = %f\n", n, h_gaussian_1D[i]);
            #endif // SHOW_TAP_VALUES
            dc_gain_1D += h_gaussian_1D[i];
        }

        // // Hack to normalize Gaussian to have DC gain = 1
        // float dc_gain_1D_hack = 0.0;
        // for (int i = 0; i < DIM; i++) {
        //     h_gaussian_1D[i] /= dc_gain_1D;
        //     dc_gain_1D_hack += h_gaussian_1D[i];
        // }

        float dc_gain = 0.0;
        for (int i = 0; i < DIM; i++) {
            // dc_gain += (dc_gain_1D_hack * h_gaussian_1D[i]);
            dc_gain += (dc_gain_1D * h_gaussian_1D[i]);
        }

        printf("h1 region of support is [-%d, %d]^2\n", H, 1);
        printf("h2 region of support is [-%d, %d]^2\n", 1, H);
        printf("DC gain of 1D Gaussian is %f\n", dc_gain_1D);
        // printf("DC gain of 1D Gaussian normalized is %f\n", dc_gain_1D_hack);
        printf("Total DC gain is %f\n", dc_gain);

        // Use cascaded separable Gaussian filters
        read_image_from_bmp(&imageIn, &input_bmp, H);
        perform_boundary_extension(&imageIn);
        apply_separable_filters(&imageIn, &imageLowPass, 
            h_gaussian_1D, h_gaussian_1D, H, H);

    }

    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageLowPass, "out_task1_lowpass.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS

    ////////////////////////////////////////////////////////
    // Task 2: Approximate gradient via finite difference //
    ////////////////////////////////////////////////////////
    
    // Currently, image borders are removed after filtering operations.
    // We need to create a new image with a border before each filtering step.
    // TODO: Can we make this more efficient?
    //
    //                +--------+                  +----------+
    // imageLowPass ->| border |-> imageGradIn -->| gradient |--> imageGrad
    //                +--------+                  +----------+

    // Side note: don't try to test this with four-pixel test images. The diff
    // of the odd symmetric border extended image will be zero...

    image imageGradIn, imageGrad;

    pixel_t hD_1[3] = {0.5, 0.0, -0.5};   // Hey, this has dc gain of 0...
    pixel_t hD_2[3] = {0.5, 0.0, -0.5};

    // Scale by alpha. This is necessary to "brighten" the result.
    // If output image seems too dark, try increase alpha (e.g. 10.0)
    for (int tap = 0; tap < 3; tap++) {
        hD_1[tap] *= alpha;
        hD_2[tap] *= alpha;
    }
    
    copy_image(&imageLowPass, &imageGradIn, 1);
    perform_boundary_extension(&imageGradIn);
    apply_separable_filters(&imageGradIn, &imageGrad, hD_1, hD_2, 1, 1);

    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageGrad, "out_task2_grad.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS

    //////////////////////////////////////////
    // Task 2 Bonus: Derivative of Gaussian //
    //////////////////////////////////////////

    // Further reading:
    // https://hannibunny.github.io/orbook/preprocessing/04gaussianDerivatives.html
    // Course notes Ch 8 Page 7

    // The separable Gaussian low-pass filter from Task 1 combined with the 
    // finite difference derivative approximation from Task 2 to perform noise
    // reduction and edge detection respectively.
    //
    // An alternative approach is to combine the Gaussian low-pass filtering 
    // and differentiation into a single "Derivative of Gaussian" (DoG) filter.

    // This section implements a DoG filter, which can be enabled using the
    // TODO: decide on a flag for dog. How about -d?

    //               +---------------+
    // imageDogIn -> | Separable DoG | -> imageDog
    //               +---------------+

    puts("\nTask 2 bonus: DoG filter");
    
    image imageDogIn, imageDog;

    const int DOG_H = 3*sigma;  // Choose region of support [-H, +H]^2
    const int DOG_DIM = (2*DOG_H+1);

    // Design separable DoG filter
    pixel_t h_dog_1D[DOG_DIM];
    float dc_gain_dog_1D = 0.0;
    for (int i = 0; i < DOG_DIM; i++) {
        const int n = i - DOG_H;
        const float scale_a = sqrt(alpha);
        const float scale_d = - n / sqrt(sigma * sigma);
        const float scale_g = 1.0 / sqrt(2 * M_PI * sigma * sigma);
        // Note: We use sqrt because filter is applied twice.
        h_dog_1D[i] = (scale_a * scale_d * scale_g) 
            * exp(-n*n / (2*sigma*sigma));
        
        #ifdef SHOW_TAP_VALUES
        printf("h[%d] = %f\n", n, h_dog_1D[i]);
        #endif // SHOW_TAP_VALUES

        dc_gain_dog_1D += h_dog_1D[i];
    }

    // Check dc gain
    float dc_gain_dog_2D = 0.0;
    for (int i = 0; i < DOG_DIM; i++) {
        dc_gain_dog_2D += (dc_gain_dog_1D * h_dog_1D[i]);
    }

    printf("h1 region of support is [-%d, %d]^2\n", DOG_H, 1);
    printf("h2 region of support is [-%d, %d]^2\n", 1, DOG_H);
    printf("DC gain of 1D DoG is %f\n", dc_gain_dog_1D);
    printf("Total DC gain is %f\n", dc_gain_dog_2D);

    copy_image(&imageIn, &imageDogIn, DOG_H);
    // apply_filter(&imageDogIn, &imageDog, h_dog_1D, DOG_H, 1);
    apply_separable_filters(&imageDogIn, &imageDog, 
        h_dog_1D, h_dog_1D, DOG_H, DOG_H);

    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageDog, "out_task2_bonus_dog.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS

    ///////////////////////
    // Task 3: Laplacian //
    ///////////////////////



    //////////////////////////////////////////
    // Task 3a Bonus: Laplacian of Gaussian //
    //////////////////////////////////////////



    // TODO: Handle -n flag.
    // TODO: export_image_as_bmp(..., outputFile);

    return EXIT_SUCCESS;
}
