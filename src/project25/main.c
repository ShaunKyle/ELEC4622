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
void create_laplacian_subimages(image *image_in, image *subimages, int D, int H);

// CLI help message (usage, description, options list)
const char CLI_HELP[] = "\
Usage: project25 [options] <file1> <file2> <D> <H>\n\
\n\
Snitches get stitches. \n\
\n\
Option      Description\n\
-o <file>   Specify output file (default: output.bmp)\n\
-l          Stitch in Laplacian domain instead of image domain\n\
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
    char * inputFile1 = NULL;    // <inputFile>
    char * inputFile2 = NULL;    // <inputFile>
    int D = 0;                  // <D>
    int H = 0;                  // <H>
    char * outputFile = NULL;           // -o <outputFile>
    bool laplacianDomainFlag = false;   // -l

    // Parse command-line arguments
    enum { O, L, NOT_OPTION } currentArg, prevArg = NOT_OPTION;
    for (int i_arg = 1; i_arg < argc; i_arg++) {
        // Check if current argument is an option flag. If so, then the next 
        // argument will contain the value of that option.
        if (strcmp(argv[i_arg], "-o") == 0) {
            currentArg = O;
        }
        else if (strcmp(argv[i_arg], "-l") == 0) {
            currentArg = L;
            laplacianDomainFlag = true;
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
        // If the first argument is not a flag, then it must be <file1>.
        // The second argument will be <file2>.
        // The third argument will be <D>.
        // The fourth argument will be <H>.
        //
        // Case 2: 
        // If the 4th-to-last argument is not an option value or a flag, then 
        // it must be the input files (final args will be <D> and <H>).
        if (inputFile1 == NULL) {
            if ((i_arg == 1) && (currentArg == NOT_OPTION)) {
                // puts("case 1");
                inputFile1 = argv[i_arg];
                inputFile2 = argv[i_arg+1];
                D = atoi(argv[i_arg+2]);
                H = atoi(argv[i_arg+3]);
            }
            else if ((i_arg == (argc - 4)) && 
            (prevArg == NOT_OPTION) && 
            (currentArg == NOT_OPTION)) {
                // puts("case 2");
                inputFile1 = argv[i_arg];
                inputFile2 = argv[i_arg+1];
                D = atoi(argv[i_arg+2]);
                H = atoi(argv[i_arg+3]);
            }
        }

        // Store state of current argument
        prevArg = currentArg;
    }

    // Print help if no input file was provided
    if ((inputFile1 == NULL) || (inputFile2 == NULL)) {
        puts(CLI_HELP);
        return EXIT_FAILURE;
    }

    // Give output file a default name if none provided
    if (outputFile == NULL) {
        outputFile = (char *) "out.bmp";
    }
    
    // Hey look, it works!
    printf("i1 = %s\n", inputFile1);
    printf("i2 = %s\n", inputFile2);
    printf("o = %s\n\n", outputFile);

    printf("D = %d\n", D);
    printf("H = %d\n", H);
    printf("L = %s\n", laplacianDomainFlag ? 
        "Laplacian domain" : "Image domain");

    // Load input bitmap file
    bmp input1_bmp, input2_bmp;
    int fileErr = load_bmp(&input1_bmp, inputFile1);
    if (fileErr != 0) {
        print_bmp_file_error(fileErr);
        return EXIT_FAILURE;
    }
    fileErr = load_bmp(&input2_bmp, inputFile2);
    if (fileErr != 0) {
        print_bmp_file_error(fileErr);
        return EXIT_FAILURE;
    }

    // Save original height for Task 3
    const int original_height = input1_bmp.rows;
    
    printf("Info about %s\n", inputFile1);
    printf("Width x Height: %dx%d px\n", input1_bmp.cols, input1_bmp.rows);
    printf("Components: %d\n", input1_bmp.num_components);
    printf("Line bytes: %d\n", input1_bmp.line_bytes);
    printf("Alignment padding bytes: %d\n\n", input1_bmp.alignment_bytes);

    printf("Info about %s\n", inputFile2);
    printf("Width x Height: %dx%d px\n", input2_bmp.cols, input2_bmp.rows);
    printf("Components: %d\n", input2_bmp.num_components);
    printf("Line bytes: %d\n", input2_bmp.line_bytes);
    printf("Alignment padding bytes: %d\n\n", input2_bmp.alignment_bytes);


    ///////////////////////////////////////
    // Task 5: Stitching in image domain //
    ///////////////////////////////////////

    // Stitching border will be located at middle of image.
    // Assuming both images have the same dimensions.

    if (!laplacianDomainFlag) {
        image image1, image2;
        read_image_from_bmp(&image1, &input1_bmp, 0);
        read_image_from_bmp(&image2, &input2_bmp, 0);

        // Allocate memory for stitched image
        image stitched;
        const int planes = image1.num_components;
        stitched.rows = image1.rows;
        stitched.cols = image1.cols;
        stitched.border = 0;
        stitched.stride = stitched.cols + 2*stitched.border;
        stitched.num_components = image1.num_components;
        stitched.handle = malloc(
            stitched.rows * stitched.cols * stitched.num_components
        );
        stitched.buf = stitched.handle;

        // Place image 1 on left
        for (int r = 0; r < image1.rows; r++) {
            // Image 1 on left
            for (int c = 0; c < (image1.cols/2); c++) {
                const int index1 = (r*image1.stride + c) * planes;
                const int indexS = (r*stitched.stride + c) * planes;

                for (int p = 0; p < planes; p++) {
                    stitched.buf[indexS+p] = image1.buf[index1+p];
                }
            }

            // // Image 1 and 2 stitching
            // for (int c = 0; c < (image1.cols/2); c++) {
            //     const int c1 = c + image1.cols/2;
            //     const int index1 = (r*image1.stride + c1) * planes;
            //     const int index2 = (r*image2.stride + c) * planes;
            //     const int sc = c + image1.cols/2;
            //     const int indexS = (r*stitched.stride + sc) * planes;

            //     for (int p = 0; p < planes; p++) {
            //         stitched.buf[indexS+p] = (image1.buf[index1+p] + image2.buf[index2+p]) / 2;
            //     }
                
            // }

            // Image 2 on right
            for (int c = 0; c < (image1.cols/2); c++) {
                const int c2 = c + image1.cols/2;
                const int index2 = (r*image2.stride + c2) * planes;
                const int sc = c + image1.cols/2;
                const int indexS = (r*stitched.stride + sc) * planes;

                for (int p = 0; p < planes; p++) {
                    stitched.buf[indexS+p] = image2.buf[index2+p];
                }
            }
        }

        // Output
        export_image_as_bmp(&stitched, outputFile);

        return EXIT_SUCCESS;
    }
    
    
    
    ///////////////////////////////////////////////////////
    // Task 2: Create Laplacian Pyramids of input images //
    ///////////////////////////////////////////////////////

    // This is a prerequisite step for Task 6.

    // Modification for Task 4: 
    // Save Laplacian detail images for reconstruction
    image image1, image2;
    read_image_from_bmp(&image1, &input1_bmp, 0);
    read_image_from_bmp(&image2, &input2_bmp, 0);

    image laplacian1Levels[D+1];
    image laplacian2Levels[D+1];
    create_laplacian_subimages(&image1, laplacian1Levels, D, H);
    create_laplacian_subimages(&image2, laplacian2Levels, D, H);

    
    // // Debug only.
    // // Output, level shift by half max magnitude to illustrate +/- values
    // for (int level = 0; level<D; level++) {
    //     perform_level_shift(&laplacian2Levels[level], 0.5);
    //     char str[20];
    //     sprintf(str, "out_x1_%d.bmp", level);
    //     export_image_as_bmp(&laplacian2Levels[level], str);
    // }
    


    ///////////////////////////////////////////////
    // Task 6: Stitch images in Laplacian domain //
    ///////////////////////////////////////////////

    // How to stitch in Laplacian domain
    // 1. Stitch together Laplacian pyramid sub-images
    // 2. Perform image reconstruction
    //
    // The following is Step 1.

    image stitchedLapLevels[D+1];

    for (int level = 0; level<D+1; level++) {
        image img1 = laplacian1Levels[level];
        image img2 = laplacian2Levels[level];

        // Allocate memory for stitchedLapLevels[level] image
        const int planes = img1.num_components;
        stitchedLapLevels[level].rows = img1.rows;
        stitchedLapLevels[level].cols = img1.cols;
        stitchedLapLevels[level].border = 0;
        stitchedLapLevels[level].stride = stitchedLapLevels[level].cols + 2*stitchedLapLevels[level].border;
        stitchedLapLevels[level].num_components = img1.num_components;
        stitchedLapLevels[level].handle = malloc(
            stitchedLapLevels[level].rows * stitchedLapLevels[level].cols * stitchedLapLevels[level].num_components
        );
        stitchedLapLevels[level].buf = stitchedLapLevels[level].handle;

        printf("planes=%d\n", stitchedLapLevels[level].num_components);
        printf("rows=%d\n", stitchedLapLevels[level].rows);

        
        // Place image 1 on left
        for (int r = 0; r < img1.rows; r++) {
            // Image 1 on left
            for (int c = 0; c < (img1.cols/2); c++) {
                const int index1 = (r*img1.stride + c) * planes;
                const int indexS = (r*stitchedLapLevels[level].stride + c) * planes;

                for (int p = 0; p < planes; p++) {
                    stitchedLapLevels[level].buf[indexS+p] = img1.buf[index1+p];
                }
            }

            // Image 2 on right
            for (int c = 0; c < (img1.cols/2); c++) {
                const int c2 = c + img1.cols/2;
                const int index2 = (r*img2.stride + c2) * planes;
                const int sc = c + img1.cols/2;
                const int indexS = (r*stitchedLapLevels[level].stride + sc) * planes;

                for (int p = 0; p < planes; p++) {
                    stitchedLapLevels[level].buf[indexS+p] = img2.buf[index2+p];
                }
            }
        }
        

    }

    // Debug only.
    // Output, level shift by half max magnitude to illustrate +/- values
    for (int level = 0; level<(D+1); level++) {
        perform_level_shift(&stitchedLapLevels[level], 0.5);
        char str[20];
        sprintf(str, "out_xs1_%d.bmp", level);
        export_image_as_bmp(&stitchedLapLevels[level], str);
    }

    
    ////////////////////////////////////////
    // Task 3: Reconstruct original image //
    ////////////////////////////////////////

    // Modification for Task 6: Use stitchedLapLevels pyramid instead of laplacianLevels
    

    // Modification for Task 4:
    // Use saved laplacian subimages, rather than reading from BMP
    
    // For lab demo: Multiply samples of one Laplacain detail image
    // perform_scaling(&laplacianLevels[0], 5.0);

    // return EXIT_SUCCESS;    // DEBUG TODO REMOVE
    puts("\nStep 3> Reconstruction");

    // Design windowed sinc interpolator
    // Confine bandwidth to (-pi/2, pi/2)^2
    // Because this filter is designed to interpolate from a downsampled 
    // image, we don't need to stretch the sinc.
    // Because the sample density is increasing.
    // See Ch 3 Pg 15 Sample density increase.
    puts("Separable interpolator g[n]");
    pixel_t g_interp[2*H+1];
    design_windowed_sinc_filter(g_interp, H, 1.0);

    // Info about original image
    printf("Original image height: %d\n\n", original_height);

    // // Load Laplacian pyramid BMP
    // bmp pyramid_bmp;
    // int fileErr2 = load_bmp(&pyramid_bmp, outputFile);
    // if (fileErr2 != 0) {
    //     print_bmp_file_error(fileErr2);
    //     return EXIT_FAILURE;
    // }
    // printf("Info about %s (Laplacian pyramid)\n", outputFile);
    // printf("Width x Height: %dx%d px\n",pyramid_bmp.cols,pyramid_bmp.rows);
    // printf("Components: %d\n", pyramid_bmp.num_components);
    // printf("Line bytes: %d\n", pyramid_bmp.line_bytes);
    // printf("Alignment padding bytes: %d\n\n", pyramid_bmp.alignment_bytes);

    // image imgPyramid;
    // read_image_from_bmp(&imgPyramid, &pyramid_bmp, 0);

    // Inverting the Laplacian Transform
    //
    // The figure below shows a segment of a Laplacian pyramid.
    //
    // | etc.      |
    // +-----+-----+
    // |     |
    // |     | y_{d}
    // |     |
    // +--+--+
    // |  | y_{d+1}
    // +--+
    //
    // To reconstruct x_d (original image x, level d of pyramid), we use a 
    // similar formula to the Laplacian pyramid.
    //
    // x_d = y_d + sum_k( x_{d+1}[k] g[n-2k] )
    //
    // Start with x_D = y_D, and keep iterating this until x_0 is reached.

    image imageX, imageXLowRes, imageY;
    int accessed_rows = 0;
    for (int level = D; level >= 0; level--) {
        printf("Level %d\n", level);
        // printf("Starting row: %d\n", accessed_rows);

        // Step 1: Obtain y_d from Laplacian pyramid
        const int height = original_height / pow(2, level);

        copy_image(&stitchedLapLevels[level], &imageY, 0);

        // // Debug: Check if y_d obtained correctly
        // char str[20];
        // sprintf(str, "out_lvl%d.bmp", level);
        // export_image_as_bmp(&imageY, str);


        // Step 2: Form original image x_d
        //
        // x_d = y_d + sum_k( x_{d+1}[k] g[n-2k] )

        if (level == D) {
            // For first iteration, x_D = y_D
            copy_image(&imageY, &imageX, 0);
        }
        else {
            // Obtain x_{d+1}
            copy_image(&imageX, &imageXLowRes, 0);

            // Intialize x_d = y_d
            copy_image(&imageY, &imageX, 0);

            // For each location k in x_{d+1}, modify x_d[2k]
            // k = [col, row]
            for (int row = 0; row < imageXLowRes.rows; row++) {
                for (int col = 0; col < imageXLowRes.cols; col++) {
                    // Pixel index for x_{d+1}[k]
                    // This acts as an "excitation source" for g_interp
                    const int planes = imageXLowRes.num_components;
                    const int index = 
                        (row * imageXLowRes.stride + col) * planes;
                    
                    // Pixel index for x_d[2k]
                    // This image is scaled up by a factor of 2
                    const int indexX = 2*(row * imageX.stride + col) * planes;

                    // Direct implementation
                    // y_d[2k] = x_d[2k] - x_{d+1}[k] * sinc * sinc
                    for (int plane = 0; plane < planes; plane++) {
                        // Obtain direct 2D PSF values
                        // Also "excite" the values of g_direct by x_{d+1}[k]
                        const int DIM = 2*H+1;
                        pixel_t g_direct[DIM*DIM];
                        for (int r = 0; r < DIM; r++) {
                            for (int c = 0; c < DIM; c++) {
                                g_direct[r*DIM + c] = g_interp[r] * g_interp[c] 
                                    * imageXLowRes.buf[index+plane];
                            }
                        }

                        // Apply excited PSF to output image
                        pixel_t *x_p = imageX.buf + (indexX+plane);
                        for (int r = 0; r < DIM; r++) {
                            for (int c = 0; c < DIM; c++) {
                                const int n1 = c - H;
                                const int n2 = r - H;
                                const int stride = imageX.stride;

                                x_p[(n2*stride + n1) * planes] -= 
                                    g_direct[r*DIM + c];
                            }
                        }

                        // Hack: x_d[2k+1] = x_{d+1}[k]
                        // This seems wrong, but whatever.
                        const int plusOneRow = 1 * planes;
                        const int plusOneCol = 1*imageX.stride*planes;
                        x_p[plusOneRow] -= imageXLowRes.buf[index+plane];
                        x_p[plusOneCol] -= imageXLowRes.buf[index+plane];
                        x_p[plusOneRow+plusOneCol] -= imageXLowRes.buf[index+plane];
                    }
                }
            }
        }

        // Debug: Check if x_d obtained correctly
        char str[20];
        sprintf(str, "out_x_%d.bmp", level);
        export_image_as_bmp(&imageX, str);

        // Ready for next row
        accessed_rows += height;
        free(imageY.handle);
        imageY.handle = NULL;
        imageY.buf = NULL;
    }

    // Step 3: Export reconstructed image, after iterating to level 0
    export_image_as_bmp(&imageX, "out_reconstruct.bmp");


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


void create_laplacian_subimages(image *image_in, image *subimages, int D, int H) {
    puts("Part 1> Laplacian Pyramid");

    // Design windowed sinc interpolator
    // Confine bandwidth to (-pi/2, pi/2)^2
    // Because this filter is designed to interpolate from a downsampled 
    // image, we don't need to stretch the sinc.
    // Because the sample density is increasing.
    // See Ch 3 Pg 15 Sample density increase.
    puts("Separable interpolator g[n]");
    pixel_t g_interp[2*H+1];
    design_windowed_sinc_filter(g_interp, H, 1.0);

    // Design windowed stretched sinc
    // Reminder: Use apply_separable_filters_2n() to apply filter
    puts("Separable stretched sinc h[n]");
    pixel_t h_sinc[2*H+1];
    design_windowed_sinc_filter(h_sinc, H, 2.0);


    // Create sub-images of pyramid
    int accessed_rows = 0;
    image imageX, imageXLowRes, imageY;
    copy_image(image_in, &imageX, H);
    // imageX = image0;
    for (int level = 0; level <= D; level++) {
        printf("Level %d\n", level);
        printf("Starting row: %d\n", accessed_rows);

        // TODO:
        // 1. Obtain x_d and x_{d+1}, where d is the current level
        // 2. Form the Laplacian detail image
        // 3. Write the Laplacian detail image to pyramid

        // Step 1: Obtain x_d and x_{d+1}
        //
        // Apply decimation filter to image x_d (every second pixel)
        // copy_image(&imageCurrent, &imagePrev, H);
        perform_boundary_extension(&imageX);
        apply_separable_filters_2n(&imageX, &imageXLowRes, 
            h_sinc, h_sinc, H, H
        );

        // Step 2: Form Laplacian detail image
        //
        // y_d[n] = x_d[n] - sum( x_{d+1}[k] g[n-2k] )
        //
        // Using input-based implementation.
        // This means that for every pixel in x_{d+1}, we add something to
        // x_d.

        // Set y_d[n] = x_d[n]
        copy_image(&imageX, &imageY, H);

        // TODO: finish.
        // For each location k in x_{d+1}, modify y_d[2k]
        // k = [col, row]
        for (int row = 0; row < imageXLowRes.rows; row++) {
            for (int col = 0; col < imageXLowRes.cols; col++) {
                // Pixel index for x_{d+1}[k]
                // This acts as an "excitation source" for g_interp
                const int planes = imageXLowRes.num_components;
                const int index = 
                    (row * imageXLowRes.stride + col) * planes;
                
                // Pixel index for y_d[2k]
                // This image is scaled up by a factor of 2
                const int indexY = 2*(row * imageY.stride + col) * planes;

                // Estimate x_d[2k] from x_{d+1}[k]
                // Note that this is a sequence, not a single value?

                // Direct implementation
                // y_d[2k] = x_d[2k] - x_{d+1}[k] * sinc * sinc
                for (int plane = 0; plane < planes; plane++) {
                    // Obtain direct 2D PSF values
                    // Also "excite" the values of g_direct by x_{d+1}[k]
                    const int DIM = 2*H+1;
                    pixel_t g_direct[DIM*DIM];
                    for (int r = 0; r < DIM; r++) {
                        for (int c = 0; c < DIM; c++) {
                            g_direct[r*DIM + c] = g_interp[r] * g_interp[c] 
                                * imageXLowRes.buf[index+plane];
                        }
                    }

                    // Apply excited PSF to output image
                    pixel_t *y_p = imageY.buf + (indexY+plane);
                    for (int r = 0; r < DIM; r++) {
                        for (int c = 0; c < DIM; c++) {
                            const int n1 = c - H;
                            const int n2 = r - H;
                            const int stride = imageY.stride;

                            y_p[(n2*stride + n1) * planes] -= 
                                g_direct[r*DIM + c];
                        }
                    }

                    // Hack: x_d[2k+1] = x_{d+1}[k]
                    // This seems wrong, but whatever.
                    const int plusOneRow = 1 * planes;
                    const int plusOneCol = 1*imageY.stride*planes;
                    y_p[plusOneRow] -= imageXLowRes.buf[index+plane];
                    y_p[plusOneCol] -= imageXLowRes.buf[index+plane];
                    y_p[plusOneRow+plusOneCol] -= imageXLowRes.buf[index+plane];
                }

            }
        }

        // Step 3: Write sub-image to pyramid
        // Modification for Task 4: Also save subimage
        copy_image(&imageY, &subimages[level], 0);
        for (int row = 0; row < imageY.rows; row++) {
            for (int col = 0; col < imageY.cols; col++) {
                // // Pixel to read from imageY
                // const int row_flipped = imageY.rows - row - 1;
                // const int indexCurrent = 
                //     (row_flipped * imageY.stride + col) 
                //     * imageY.num_components;
                
                // // Calculate pyramid coordinates and array index
                // // Flip n2 because of BMP coordinate system
                // const int n1 = col;
                // const int n2 = total_height - accessed_rows - 1;
                // const int indexPyramid = 
                //     (n2 * imgLaplacianPyramid.stride + n1) 
                //     * imgLaplacianPyramid.num_components;
                
                // // Write pixel from imageY to pyramid
                // for (int plane = 0; plane < imageY.num_components; plane++) {
                //     imgLaplacianPyramid.buf[indexPyramid+plane] = 
                //         imageY.buf[indexCurrent+plane];
                // }
            }

            // Ready to write next row
            // printf("%d\n", accessed_rows);
            accessed_rows++;
        }

        // Ready for next level of pyramid
        copy_image(&imageXLowRes, &imageX, H);
    }


}
