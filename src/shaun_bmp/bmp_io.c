#include <stdio.h>
#include <stdint.h>
#include <string.h> // memset
#include <assert.h>

#include "bmp_io.h"

// Bitmap header memory offsets
const int size_offset = 0x02;
const int start_address_offset = 0x0A;

// Convert value from little-endian to whatever machine endian is.
// Assumes that word size is 32-bits.
static void from_little_endian(int32_t *words, int num_words) {
    // Check machine endianness (in C++20, we can just use std::endian)
    int test = 1;
    unsigned char *first_byte = (unsigned char *) &test;
    if (*first_byte == 1) {
        // Machine uses little-endian representation, no change needed.
        // This is most likely the case on personal computers.
        return;
    }

    // Machine uses big-endian representation, change order of bytes.
    int32_t tmp;
    for (; num_words--; words++) {
        tmp = *words;
        *words = ((tmp >> 24) & 0x000000FF) +
                 ((tmp >> 8)  & 0x0000FF00) +
                 ((tmp << 8)  & 0x00FF0000) +
                 ((tmp << 24) & 0xFF000000);
    }
}

// Convert value from machine-endian to little-endian
// Assumes that word size is 32-bits.
inline static void to_little_endian(int32_t *words, int num_words) {
    from_little_endian(words, num_words);
}

//! \brief Load a bitmap file for read/write operations.
//!
//! Opens the file, extracts relevant header info, then populates a bmp struct.
//!
//! @param bmp_info Bitmap struct to populate with information about file
//! @param fname    Name/path of bitmap file to load
//!
//! @return Bitmap IO_ERR
int load_bmp(bmp *bmp_info, const char *fname) {
    // Open binary file for read operations. Error if file cannot be opened.
    FILE *in = fopen(fname,"rb");
    if (in == NULL) {
        return IO_ERR_NO_FILE;
    }

    ///////////////////////////////////////
    // File header (info about the file) //
    ///////////////////////////////////////

    // Read bitmap file header (first 14 bytes)
    uint8_t magic[14];
    fread(magic, 1, 14, in);

    // Check signature
    if (magic[0] != 'B' || magic[1] != 'M') {
        fclose(in);
        return IO_ERR_FILE_HEADER;
    }

    // Extract file size and starting address of image data (32-bit int)
    //
    // Note that these values are stored in little-endian format. We have two 
    // options for reading these values:
    //
    // Option 1: Read each byte individually to reconstruct the 32-bit integer 
    //           with the byte order reversed.
    //
    // int32_t size = 0;
    // for (int byte = 0; byte < 4; byte++) {
    //     const int memory_offset = 2;
    //     size |= (magic[byte+memory_offset] << 8*byte);
    // }
    //
    // Option 2: Copy entire 32-bit memory contents, only reverse the order of 
    //           bytes if the machine memory architecture uses big-endian.

    int32_t *size_ptr = (int32_t *) &magic[size_offset];
    int32_t size = *size_ptr;
    from_little_endian(&size, 1);   // This won't do anything on a machine that 
                                    // uses little-endian format already.

    int32_t *start_address_ptr = (int32_t *) &magic[start_address_offset];
    int32_t start_address = *start_address_ptr;
    from_little_endian(&start_address, 1);

    // Printf debugging...
    // puts("\nFile header info");
    // printf("size = %d bytes\n", size);
    // printf("image start address (offset) = %d\n", start_address);

    ///////////////////////////////////////
    // DIB header (info about the image) //
    ///////////////////////////////////////

    // Windows BITMAPINFOHEADER header

    // Read the bitmap information header  (DIB header)
    // 10 items * 4 bytes each
    bmp_header header;
    size_t dib_header_bytes = fread(&header, 1, 40, in);
    if (dib_header_bytes != 40) {
        fclose(in);
        return IO_ERR_FILE_TRUNC;
    }
    from_little_endian((int32_t *) &header, 10);

    // Check reported size of DIB header is 40 B (10 items * 4 bytes each)
    if (header.size != 40) {
        fclose(in);
        return IO_ERR_FILE_HEADER;
    }

    // Check no. of planes is 1
    if (header.planes != 1) {
        fclose(in);
        return IO_ERR_FILE_HEADER;
    }

    // Figure out no. of color components based on the color depth
    // We only deal with 24-bit RGB and 8-bit greyscale currently.
    int num_components;
    if (header.bit_count == 24) {
        num_components = 3;
    }
    else if (header.bit_count == 8) {
        num_components = 1;
    }
    else {
        fclose(in);
        return IO_ERR_UNSUPPORTED;
    }

    // Check if compression type is supported by this library
    // Currently, we only work with uncompressed files.
    if (header.compression != 0) {
        fclose(in);
        return IO_ERR_UNSUPPORTED;
    }

    // TODO: Color palette

    // Printf debugging...
    // puts("\nDIB header info");
    // printf("Size of header (should be 40) = %d bytes\n", header.size);
    // printf("Size of image (could be 0 if BI_RGB) = %d\n", header.image_size);
    // printf("No. of components = %d\n", num_components);
    // printf("Width and Height = %dx%d px\n", header.width, header.height);
    // printf("Compression method = %d\n", header.compression);

    // Populate bmp_info fields
    bmp_info->rows = header.height;
    bmp_info->cols = header.width;
    bmp_info->num_components = num_components;
    bmp_info->line_bytes = bmp_info->num_components * bmp_info->cols;
    bmp_info->alignment_bytes = 
        (4 - bmp_info->line_bytes) & 3; // Pad to a multiple of 4 bytes
    bmp_info->num_unaccessed_rows = bmp_info->rows;
    bmp_info->file = in;

    return 0;
}

//! \brief Create a new bitmap file to write image data to.
//!
//! Opens a new file, writes a bitmap header, and populates a bmp struct to 
//! access the file for future write operations.
//!
//! @return Bitmap IO_ERR
int create_bmp(bmp *bmp_info, const char *fname, int width, int height, 
               int num_components) {
    // Open new file to write bitmap data to.
    FILE *out = fopen(fname, "wb");
    if (out == NULL) {
        return IO_ERR_NO_FILE;
    }

    // Calculate theoretical header size, file size
    int header_size_bytes = 14 + sizeof(bmp_header); // File header = 14 B
    assert(header_size_bytes == 54);                 // DIB = 10 * 4 B = 40 B
    if (num_components == 1) {
        header_size_bytes += 1024;  // Need a color lookup table
    } else if (num_components != 3) {
        fclose(out);
        return IO_ERR_UNSUPPORTED;  // Only support 1 or 3 components
    }
    int file_size_bytes = header_size_bytes + 
        (bmp_info->line_bytes + bmp_info->alignment_bytes) * bmp_info->rows;

    // Write bitmap file header
    uint8_t magic[14];
    magic[0] = 'B'; magic[1] = 'M'; // File signature "BM"
    for (int n = 0; n < 4; n++) {   // File size
        magic[size_offset + n] = (uint8_t) (file_size_bytes >> 8 * n);
    }
    magic[6] = magic[7] = magic[8] = magic[9] = 0;
    for (int n = 0; n < 4; n++) {   // Image data start address (after header)
        magic[start_address_offset + n] = (uint8_t) (header_size_bytes >> 8 * n);
    }
    fwrite(magic, 1, 14, out);

    // Write bitmap DIB header
    bmp_header header;
    header.size = 40;
    header.width = width;
    header.height = height;
    header.planes = 1;
    header.bit_count = (num_components==1)?8:24;
    header.compression = 0; // BI_RGB (no compression)
    header.image_size = 0;  // Using BI_RGB, can be 0
    header.xpels_per_metre = header.ypels_per_metre = 0;
    header.num_colours_used = header.num_colours_important = 0;
    to_little_endian((int32_t *) &header, 10);
    fwrite(&header, 1, 40, out);

    // Write color LUT if required
    if (num_components == 1) {
        for (int n = 0; n < 256; n++) {
            fputc(n, out);
            fputc(n, out);
            fputc(n, out);
            fputc(0, out);
        }
    }

    // Reset bmp_info struct, populate bmp_info
    memset(bmp_info, 0, sizeof(bmp));
    bmp_info->rows = bmp_info->num_unaccessed_rows = height;
    bmp_info->cols = width;
    bmp_info->num_components = num_components;
    bmp_info->line_bytes = num_components * width;
    bmp_info->alignment_bytes = (4 - bmp_info->line_bytes) & 3;
    bmp_info->file = out;

    return 0;
}

//! \brief Close a bitmap file and reset the associated bmp struct.
//!
//! @param bmp_info Bitmap struct with file handle to close
//!
//! @return Bitmap IO_ERR
void close_bmp(bmp *bmp_info) {
    // Check if there is a file to close
    if (bmp_info->file != NULL) {
        fclose(bmp_info->file);
    }

    // Reset bmp_info struct
    memset(bmp_info, 0, sizeof(bmp));
    bmp_info->file = NULL;
}

//! \brief Read the next line in the bitmap image
//!
//! @param bmp_info Bitmap struct with file handle to read from
//! @param line     Buffer to read line data from
//!
//! @return Bitmap IO_ERR
int read_bmp_line(bmp *bmp_info, uint8_t *line) {
    // Precondition: There is a line available to read.
    if ((bmp_info->file == NULL) || (bmp_info->num_unaccessed_rows <= 0)) {
        return IO_ERR_FILE_NOT_OPEN;
    }

    // Read next line from file
    size_t line_bytes_read = 
        fread(line, 1, (size_t) bmp_info->line_bytes, bmp_info->file);
    if (bmp_info->alignment_bytes > 0) {
        uint8_t buf[3]; // Padding can be 0 to 3 bytes. (add to multiple of 4)
        fread(buf, 1, (size_t) bmp_info->alignment_bytes, bmp_info->file);
    }
    bmp_info->num_unaccessed_rows--;

    // Post condition: Bytes read matches line_bytes
    if (line_bytes_read != (size_t) bmp_info->line_bytes) {
        return IO_ERR_FILE_TRUNC;
    }

    return 0;
}

//! \brief Write the next line in the bitmap image
//!
//! @param bmp_info Bitmap struct with file handle to write to
//! @param line     Buffer to write line data to
//!
//! @return Bitmap IO_ERR
int write_bmp_line(bmp *bmp_info, uint8_t *line) {
    // Precondition: There is a line available to write to.
    if ((bmp_info->file == NULL) || (bmp_info->num_unaccessed_rows <= 0)) {
        return IO_ERR_FILE_NOT_OPEN;
    }

    // Write line to next row in file
    size_t line_bytes_written = 
        fwrite(line, 1, (size_t) bmp_info->line_bytes, bmp_info->file);
    if (bmp_info->alignment_bytes > 0) {
        uint8_t buf[3] = {0,0,0}; // Padding can be 0 to 3 bytes
        fwrite(buf, 1, (size_t) bmp_info->alignment_bytes, bmp_info->file);
    }
    bmp_info->num_unaccessed_rows--;

    // Postcondition: Bytes written matches line_bytes
    if (line_bytes_written != (size_t) bmp_info->line_bytes) {
        return IO_ERR_FILE_TRUNC;
    }

    return 0;
}

//! \brief Print information about error code to stderr
void print_bmp_file_error(int fileErr) {
    if (fileErr == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
    else if (fileErr == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
    else if (fileErr == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
    else if (fileErr == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
    else if (fileErr == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
}
