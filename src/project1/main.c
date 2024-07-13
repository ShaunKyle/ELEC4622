#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#define _USE_MATH_DEFINES
#include <math.h>   // fabs, exp, M_PI

#include "../shaun_bmp/bmp_io.h"
#include "../shaun_bmp/image.h"

// Debug flags
#define EXPORT_INTERMEDIATE_STEPS
// #define SHOW_TAP_VALUES

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: project1 [options] <file> <sigma> <alpha> <beta>\n\
\n\
Performs low pass filtering and differentiation on image. \n\
Sigma is std dev between 1.0 to 100.0\n\
Alpha is a real-valued positive scaling factor.\n\
Beta is a real-valued positive scaling factor.\n\
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
    float beta = 0.0;           // <beta>
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
        // The fourth argument will be beta.
        //
        // Case 2: 
        // If the 4th-to-last argument is not an option value or a flag, then 
        // it must be the input file (final args will be sigma, alpha, beta).
        if (inputFile == NULL) {
            if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
                // puts("case 1");
                inputFile = argv[i_arg];
                sigma = atof(argv[i_arg+1]);
                alpha = atof(argv[i_arg+2]);
                beta = atof(argv[i_arg+3]);
            }
            else if ((i_arg == (argc - 4)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
                // puts("case 2");
                inputFile = argv[i_arg];
                sigma = atof(argv[i_arg+1]);
                alpha = atof(argv[i_arg+2]);
                beta = atof(argv[i_arg+3]);
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
    printf("beta = %f\n", beta);
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

        // Task 1 bonus
        // apply_optimized_moving_average_filter(&imageIn, &imageLowPass, H, H);
        
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

    const int DIFF_H = 1;
    // const int DIFF_DIM = (2*DIFF_H+1);
    #define DIFF_DIM    3

    pixel_t hD_1[DIFF_DIM] = {0.5, 0.0, -0.5}; // Hey, this has dc gain of 0...
    pixel_t hD_2[DIFF_DIM] = {0.5, 0.0, -0.5};

    // // Scale by alpha. This is necessary to "brighten" the result.
    // // If output image seems too dark, try increase alpha (e.g. 10.0)
    // for (int tap = 0; tap < DIFF_DIM; tap++) {
    //     // Scale final image by alpha by scaling filter tap values.
    //     // alpha * 1 = alpha
    //     hD_1[tap] *= alpha;
    //     hD_2[tap] *= 1;
    // }
    
    copy_image(&imageLowPass, &imageGradIn, DIFF_H);
    perform_boundary_extension(&imageGradIn);
    apply_separable_filters(&imageGradIn, &imageGrad, hD_1, hD_2, DIFF_H, DIFF_H);

    // Scale final image by alpha
    perform_scaling(&imageGrad, alpha);

    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageGrad, "out_task2_grad.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS


    //////////////////////////////////////////
    // Task 2 Bonus: Derivative of Gaussian //
    //////////////////////////////////////////

    // Further reading:
    // https://hannibunny.github.io/orbook/preprocessing/04gaussianDerivatives.html
    // Course notes Ch 8 Page 7
    // Ch 8 pg 8 Section 2.1.5 Efficient Implementation of DOG filters

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
        // const float scale_a = sqrt(alpha);
        const float scale_d = - n / (sigma * sigma);
        const float scale_g = 1.0 / sqrt(2 * M_PI * sigma * sigma);
        h_dog_1D[i] = (1 * scale_d * scale_g) 
            * exp(-n*n / (2*sigma*sigma));
        
        #ifdef SHOW_TAP_VALUES
        printf("h[%d] = %f\n", n, h_dog_1D[i]);
        #endif // SHOW_TAP_VALUES

        dc_gain_dog_1D += h_dog_1D[i];
    }

    // TODO: Try Difference of Gaussian. More efficient?
    // See Ch 8 Section 2.1.5
    // // 1D Gaussian
    // pixel_t h_gaussian_1D[DOG_DIM];
    // for (int i = 0; i < DOG_DIM; i++) {
    //     const int n = i - DOG_H;
    //     const float scale = 1.0 / (2 * M_PI * sigma * sigma);
    //     h_gaussian_1D[i] = sqrt(scale) * exp(-n*n / (2*sigma*sigma));
    // }

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
    perform_boundary_extension(&imageDogIn);
    // apply_filter(&imageDogIn, &imageDog, h_dog_1D, DOG_H, 1);
    apply_separable_filters(&imageDogIn, &imageDog, 
        h_dog_1D, h_dog_1D, DOG_H, DOG_H);
    
    // Scale final image intensity by alpha
    perform_scaling(&imageDog, alpha);

    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageDog, "out_task2_bonus_dog.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS


    ///////////////////////////////////
    // Task 3: Approximate Laplacian //
    ///////////////////////////////////

    // Step 1: High frequency noise reduction (from Task 1)
    //
    //             +----------+
    // imageIn --> | Low pass | --> imageLowPass (copy to imageLapIn)
    //             +----------+
    //
    // Step 2: Laplacian (by applying finite difference twice)
    //
    //                +------+                    +------+
    // imageLapIn --> | diff | -> imageLapDiff -> | diff | --> imageLap
    //                +------+                    +------+

    // Desired output image is scaled by alpha and level shifted by 1/2 maximum 
    // intensity (128 in byte value).
    //
    // y[n] = alpha grad^2 x[n] + 128

    image imageLapIn, imageLapDiff, imageLapDiffCopy, imageLap;

    // Design finite difference filter for Laplacian
    pixel_t hD_lap_1[DIFF_DIM] = {0.5, 0.0, -0.5};
    pixel_t hD_lap_2[DIFF_DIM] = {0.5, 0.0, -0.5};

    // // Scale final image by alpha
    // for (int tap = 0; tap < DIFF_DIM; tap++) {
    //     // TODO: Should this be alpha/2 or sqrt(alpha)?
    //     hD_lap_1[tap] *= alpha/2;
    //     hD_lap_2[tap] *= 1;
    // }

    copy_image(&imageLowPass, &imageLapIn, DIFF_H);
    perform_boundary_extension(&imageLapIn);
    apply_separable_filters(&imageLapIn, &imageLapDiff, 
        hD_lap_1, hD_lap_2, DIFF_H, DIFF_H);
    copy_image(&imageLapDiff, &imageLapDiffCopy, DIFF_H);
    perform_boundary_extension(&imageLapDiffCopy);
    apply_separable_filters(&imageLapDiffCopy, &imageLap, 
        hD_lap_1, hD_lap_2, DIFF_H, DIFF_H);
    
    // Scale final image by alpha
    perform_scaling(&imageLap, alpha);

    // Level shift final image by 1/2 maximum intensity
    perform_level_shift(&imageLap, 0.5);
    
    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageLapDiff, "out_Task3_Gradient.bmp");
    export_image_as_bmp(&imageLap, "out_Task3_Laplacian.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS


    //////////////////////////////////////////
    // Task 3a Bonus: Laplacian of Gaussian //
    //////////////////////////////////////////

    // Pretty much the same as what we did for Derivative of Gaussian, but take 
    // the second derivative of the Gaussian.

    //                   +---------------+
    // imageIn (copy) -> | Separable LoG | -> imageLog
    //                   +---------------+

    puts("\nTask 3a bonus: LoG filter");

    image imageLogIn, imageLog;

    const int LOG_H = 3*sigma;  // Choose region of support [-H, +H]^2
    const int LOG_DIM = (2*LOG_H+1);

    // Design separable LoG filter
    pixel_t h_log_1D[LOG_DIM];
    for (int i = 0; i < LOG_DIM; i++) {
        const int n = i - LOG_H;
        const float scale_dd = (n*n) / (sigma * sigma) - 1.0;
        const float scale_g = 1.0 / sqrt(2 * M_PI * sigma * sigma);
        h_log_1D[i] = (scale_dd * scale_g) * exp(-n*n / (2*sigma*sigma));

        #ifdef SHOW_TAP_VALUES
        printf("h[%d] = %f\n", n, h_log_1D[i]);
        #endif // SHOW_TAP_VALUES
    }

    printf("h1 region of support is [-%d, %d]^2\n", LOG_H, 1);
    printf("h2 region of support is [-%d, %d]^2\n", 1, LOG_H);

    copy_image(&imageIn, &imageLogIn, LOG_H);
    perform_boundary_extension(&imageLogIn);
    apply_separable_filters(&imageLogIn, &imageLog, 
        h_log_1D, h_log_1D, LOG_H, LOG_H);
    
    // // Design direct LoG filter
    // pixel_t h_log_direct[LOG_DIM*LOG_DIM];
    // for (int r = 0; r < LOG_DIM; r++) {
    //     for (int c = 0; c < LOG_DIM; c++) {
    //         const int n1 = c - LOG_H;
    //         const int n2 = r - LOG_H;
    //         const int i = c + r*LOG_DIM;
            
    //         const float term1 = 1 / (M_PI * pow(sigma, 4));
    //         const float term2 = (n1*n1 + n2*n2) / (2*sigma*sigma) - 1;
    //         h_log_direct[i] = (term1 * term2) * 
    //             exp(-(n1*n1 + n2*n2) / (2*sigma*sigma));
            
    //         printf("h[%d, %d] = %f\n", n1, n2, h_log_direct[i]);
    //     }
    // }

    // printf("h1 region of support is [-%d, %d]^2\n", LOG_H, 1);
    // printf("h2 region of support is [-%d, %d]^2\n", 1, LOG_H);

    // copy_image(&imageIn, &imageLogIn, LOG_H);
    // perform_boundary_extension(&imageLogIn);
    // apply_filter(&imageLogIn, &imageLog, h_log_direct, LOG_H, LOG_H);


    
    // Scale final image by alpha
    perform_scaling(&imageLog, alpha);

    // Level shift final image by 1/2 maximum intensity
    perform_level_shift(&imageLog, 0.5);

    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageLog, "out_task3a_bonus_log.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS


    /////////////////////////////////////////////////////////
    // Task 3b Bonus: Different operations per color plane //
    /////////////////////////////////////////////////////////

    // Blue plane:  Original image
    // Green plane: LoG scaled by beta 
    // Red plane:   DoG scaled by alpha (reuse Task 2 bonus)

    // For 3-component images, only process the green color plane

    puts("\nTask 3b: Combined");

    image imageB, imageG, imageR, imageRGB;
    image imageG_in;

    if (imageIn.num_components == 3) {
        puts("Processing green plane only");
        extract_component(&imageB, &imageIn, 1);
        copy_image(&imageB, &imageG_in, LOG_H);
        extract_component(&imageR, &imageDog, 1);
    } else {
        puts("Single plane");
        copy_image(&imageIn, &imageB, 0);
        copy_image(&imageIn, &imageG_in, LOG_H);
        copy_image(&imageDog, &imageR, 0);
    }

    // LoG scaled by beta and level shifted by 1/2 max intensity
    perform_boundary_extension(&imageG_in);
    apply_separable_filters(&imageG_in, &imageG, 
        h_log_1D, h_log_1D, LOG_H, LOG_H);
    perform_scaling(&imageG, beta);
    perform_level_shift(&imageG, 0.5);

    // Output
    #ifdef EXPORT_INTERMEDIATE_STEPS
    export_image_as_bmp(&imageIn, "out_task3b_in.bmp");
    export_image_as_bmp(&imageB, "out_task3b_blue_original.bmp");
    export_image_as_bmp(&imageG, "out_task3b_green_LoG.bmp");
    export_image_as_bmp(&imageR, "out_task3b_red_DoG.bmp");
    #endif // EXPORT_INTERMEDIATE_STEPS

    hacky_combine_planes_into_RGB(&imageRGB, &imageR, &imageG, &imageB);
    export_image_as_bmp(&imageRGB, outputFile);


    // TODO: Handle -n flag.

    return EXIT_SUCCESS;
}
