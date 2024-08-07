#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "../shaun_bmp/bmp_io.h"
#include "../shaun_bmp/image.h"
#include "../shaun_bmp/motion.h"

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
    
    printf("Info about %s\n", targetFile);
    printf("Width x Height: %dx%d px\n\n", target_bmp.cols, target_bmp.rows);
    // printf("Components: %d\n", input_bmp.num_components);
    // printf("Line bytes: %d\n", input_bmp.line_bytes);
    // printf("Alignment padding bytes: %d\n\n", input_bmp.alignment_bytes);

    //////////////////////////////////////
    // Task 1: Calculate motion vectors //
    //////////////////////////////////////

    image source, target, output;
    image output2;

    read_image_from_bmp(&source, &source_bmp, 0);
    read_image_from_bmp(&target, &target_bmp, 0);

    // SAD = sum of absolute differences
    // MSE = mean square error

    // Estimate motion vector for each block in target frame
    const int N_BLOCKS = (target.rows / B) * (target.cols / B);
    printf("No. of blocks: %d\n\n", N_BLOCKS);
    mvector_t vec[N_BLOCKS];
    mvector2_t vec2[N_BLOCKS];

    double mse[N_BLOCKS];
    double mse2[N_BLOCKS];
    for (int index = 0, r = 0; r < (target.rows - B); r += B) {
        for (int c = 0; c < (target.cols - B); c += B) {
            vec[index] = estimate_motion_block(&source, &target, r, c, B, S, 
                (minimizeMSEFlag ? MSE : MAD), &mse[index]
            );
            // printf("vec[%d] = (%d, %d)\n", index, vec[index].x, vec[index].y);

            // Task 3: Half-pel grid
            // vec2[index] = estimate_motion_block_bilinear(&source, &target, r, c, B, S, 
            //     (minimizeMSEFlag ? MSE : MAD), 4, &mse2[index]
            // );

            // // Task 4: Quarter-pel grid
            // vec2[index] = estimate_motion_block_bilinear(&source, &target, r, c, B, S, 
            //     (minimizeMSEFlag ? MSE : MAD), 4
            // );

            // // Task 5: Quarter-pel grid with windowed shifted sinc
            // vec2[index] = estimate_motion_block_sinc(&source, &target, r, c, B, S, 
            //     (minimizeMSEFlag ? MSE : MAD)
            // );

            // Task 6: Half-pel telescopic search
            vec2[index] = estimate_motion_block_telescopic(&source, &target, r, c, B, S, 
                (minimizeMSEFlag ? MSE : MAD), &mse2[index]
            );

            // Information about block motion
            // printf("Block %d\nStarting pos: [%d, %d]\n", index, c, r);
            // printf("vec(%d, %d)\n\n", vec[index].x, vec[index].y);

            // Next block
            index++;
        }
    }

    double sum_mse = 0;
    for (int n = 0; n < N_BLOCKS; n++) {
        if (mse[n] < 0)
            continue;
        if (mse[n] > 10000.0)
            continue;
        sum_mse += mse[n];
        // printf("%f\n", mse[n]);
    }
    double sum_mse2 = 0;
    for (int n = 0; n < N_BLOCKS; n++) {
        if (mse2[n] < 0)
            continue;
        if (mse2[n] > 10000.0)
            continue;
        sum_mse2 += mse2[n];
        // printf("%f\n", mse[n]);
    }
    
    printf("\nTotal image MSE: %f\n\n", sum_mse);
    printf("\nTotal bilinear2 image MSE: %f\n\n", sum_mse2);

    // Perform motion compensation on target frame
    // copy_image(&target, &target2, 0);

    image blank, blank2;
    copy_image(&target, &blank, 0);
    perform_scaling(&blank, 0.0);
    copy_image(&blank, &blank2, 0);
    for (int index = 0, r = 0; r < (target.rows - B); r += B) {
        for (int c = 0; c < (target.cols - B); c += B) {
            compensate_motion_block(&source, &blank, vec[index], r, c, B);
            compensate_motion_block_fractional(&source, &blank2, vec2[index], 
            r, c, B);
            // printf("vec[%d] = (%d, %d)\n", index, vec[index].x, vec[index].y);

            // Next block
            index++;
        }
    }

    // Draw motion vectors describing motion of each target frame block
    mono_to_RGB(&blank, &output);
    perform_scaling(&output, 0.5);
    perform_level_shift(&output, 0.5);

    mono_to_RGB(&blank2, &output2);
    perform_scaling(&output2, 0.5);
    perform_level_shift(&output2, 0.5);

    for (int index = 0, r = 0; r < (target.rows - B); r += B) {
        for (int c = 0; c < (target.cols - B); c += B) {
            // Row and col for centre of block, and end of motion vector
            const int r_centre = r + B/2;
            const int c_centre = c + B/2;
            const int r_end = r_centre + vec[index].y;
            const int c_end = c_centre + vec[index].x;

            const int r_end_fractional = r_centre + vec2[index].y/vec2[index].precision;
            const int c_end_fractional = c_centre + vec2[index].x/vec2[index].precision;
            // printf("Start: [%d, %d]\n", c_centre, r_centre);
            // printf("Vec:   (%d, %d)\n", vec[index].x, vec[index].y);
            // printf("End:   [%d, %d]\n\n", c_end, r_end);

            // Memory buf index for centre of block
            const int centre = ((r_centre * output.stride) + c_centre) * 3;
            const int end = ((r_end * output.stride) + c_end) * 3;

            const int end_fractional = ((r_end_fractional * output.stride) + c_end_fractional) * 3;

            // Draw motion vector in purple
            const int v1 = vec[index].x;
            const int v2 = vec[index].y;
            const int c1 = c_centre;
            const int c2 = r_centre;

            const int v1_fractional = vec2[index].x/vec2->precision;
            const int v2_fractional = vec2[index].y/vec2->precision;

            if (c_centre > c_end) {
                for (int n1 = c_end; n1 <= c_centre; n1++) {
                    int n2 = c2 + ((n1 - c1) * v2 / v1);
                    int i = ((n2 * output.stride) + n1) * 3;
                    output.buf[i+1] = 0;
                }
            }
            else {
                for (int n1 = c_centre; n1 <= c_end; n1++) {
                    int n2 = c2 + ((n1 - c1) * v2 / v1);
                    int i = ((n2 * output.stride) + n1) * 3;
                    output.buf[i+1] = 0;
                }
            }

            if (r_centre > r_end) {
                for (int n2 = r_end; n2 <= r_centre; n2++) {
                    int n1 = c1 + ((n2 - c2) * v1 / v2);
                    int i = ((n2 * output.stride) + n1) * 3;
                    output.buf[i+1] = 0;
                }
            }
            else {
                for (int n2 = r_centre; n2 <= r_end; n2++) {
                    int n1 = c1 + ((n2 - c2) * v1 / v2);
                    int i = ((n2 * output.stride) + n1) * 3;
                    output.buf[i+1] = 0;
                }
            }


            // Fractional line
            if (c_centre > c_end_fractional) {
                for (int n1 = c_end_fractional; n1 <= c_centre; n1++) {
                    int n2_fractional = c2 + ((n1 - c1) * v2_fractional / v1_fractional);
                    int i_fractional = ((n2_fractional * output.stride) + n1) * 3;
                    output2.buf[i_fractional+1] = 0;
                }
            }
            else {
                for (int n1 = c_centre; n1 <= c_end_fractional; n1++) {
                    int n2_fractional = c2 + ((n1 - c1) * v2_fractional / v1_fractional);
                    int i_fractional = ((n2_fractional * output.stride) + n1) * 3;
                    output2.buf[i_fractional+1] = 0;
                }
            }

            if (r_centre > r_end_fractional) {
                for (int n2 = r_end_fractional; n2 <= r_centre; n2++) {
                    int n1_fractional = c1 + ((n2 - c2) * v1_fractional / v2_fractional);
                    int i_fractional = ((n2 * output.stride) + n1_fractional) * 3;
                    output2.buf[i_fractional+1] = 0;
                }
            }
            else {
                for (int n2 = r_centre; n2 <= r_end_fractional; n2++) {
                    int n1_fractional = c1 + ((n2 - c2) * v1_fractional / v2_fractional);
                    int i_fractional = ((n2 * output.stride) + n1_fractional) * 3;
                    output2.buf[i_fractional+1] = 0;
                }
            }



            // Put red dot in centre of block
            output.buf[centre] = 0;
            output.buf[centre+1] = 0;
            output.buf[centre+2] = 1;

            output2.buf[centre] = 0;
            output2.buf[centre+1] = 0;
            output2.buf[centre+2] = 1;

            // Put green dot at end of motion vector
            output.buf[end] = 0;
            output.buf[end+1] = 1;
            output.buf[end+2] = 0;

            output2.buf[end_fractional] = 0;
            output2.buf[end_fractional+1] = 1;
            output2.buf[end_fractional+2] = 0;

            // Next block
            index++;
        }
    }

    export_image_as_bmp(&output, outputFile);
    export_image_as_bmp(&output2, "out_bilinear.bmp");

    return EXIT_SUCCESS;
}
