//! File I/O operations for bitmap files.
//!
//! Only has partial support for bitmap files (e.g. no compression)

#ifndef SHAUN_BMP_IO_H
#define SHAUN_BMP_IO_H

#include <stdint.h>
#include <stdio.h>

////////////////////////////
// Data type declarations //
////////////////////////////

typedef struct bmp_header bmp_header;
typedef struct bmp bmp;

/////////////////
// Error codes //
/////////////////

#define IO_ERR_NO_FILE          ((int) -1) // If file not found
#define IO_ERR_FILE_HEADER      ((int) -2) // If header has an error
#define IO_ERR_FILE_TRUNC       ((int) -3) // If file ends unexpectely
#define IO_ERR_UNSUPPORTED      ((int) -4) // Exception code if file uses an
                                           // unsupported format.
#define IO_ERR_FILE_NOT_OPEN    ((int) -5) // If trying to read/write a file
                                           // which is not open, or has come
                                           // to the end.

////////////////
// Structures //
////////////////

// Device independent bitmap (DIB) header
//
// Used when extracting image information from a bitmap file. 
// https://learn.microsoft.com/en-us/windows/win32/api/wingdi/ns-wingdi-bitmapinfoheader
struct bmp_header {
    uint32_t size; // Size of this structure: must be 40
    int32_t width; // Image width
    int32_t height; // Image height; -ve means top to bottom.
    // uint32_t planes_bits; // Planes in 16 LSB's (must be 1); bits in 16 MSB's
    uint16_t planes; // Number of planes (must be 1)
    uint16_t bit_count; // Number of bits per pixel

    uint32_t compression; // Only accept 0 here (uncompressed RGB data)
    uint32_t image_size; // Can be 0
    int32_t xpels_per_metre; // We ignore these
    int32_t ypels_per_metre; // We ignore these
    uint32_t num_colours_used; // Entries in colour table (0 means use default)
    uint32_t num_colours_important; // 0 means all colours are important.
};

// Bitmap image
//
// Used to store information about an open bitmap file, and access the file for 
// read/write operations.
//
// Currently, read/write operations are done sequentially row by row.
struct bmp {
    // Info
    int num_components, rows, cols;
    int line_bytes; // Number of bytes in each line, not including padding
    int alignment_bytes; // Number of 0's at end of each line.

    // State (of sequential I/O operations)
    int num_unaccessed_rows;

    // File handle
    FILE * file;
};

////////////////
// Public API //
////////////////

// Bitmap file I/O
int load_bmp(bmp *bmp_info, const char *fname);
int create_bmp(bmp *bmp_info, const char *fname, int width, int height, 
               int num_components);
void close_bmp(bmp *bmp_info);
int read_bmp_line(bmp *bmp_info, uint8_t *line);
int write_bmp_line(bmp *bmp_info, uint8_t *line);

// Error handling
void print_bmp_file_error(int fileErr);

#endif // SHAUN_BMP_IO_H
