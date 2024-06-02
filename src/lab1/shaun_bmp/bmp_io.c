#include <stdio.h>
#include <stdint.h>
#include <sys/_types/_int32_t.h>

#include "bmp_io.h"

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

// Just read the header(s) of a .bmp file.
// TODO: replace with open_bmp(), which will open fname and return an "image"
//       data structure which can be used to read image data line by line.
int read_header(const char *fname) {
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

    const int size_offset = 0x02;
    int32_t *size_ptr = (int32_t *) &magic[size_offset];
    int32_t size = *size_ptr;
    from_little_endian(&size, 1);   // This won't do anything on a machine that 
                                    // uses little-endian format already.

    const int start_address_offset = 0x0A;
    int32_t *start_address_ptr = (int32_t *) &magic[start_address_offset];
    int32_t start_address = *start_address_ptr;
    from_little_endian(&start_address, 1);

    puts("\nFile header info");
    printf("size = %d bytes\n", size);
    printf("image start address (offset) = %d\n", start_address);

    ///////////////////////////////////////
    // DIB header (info about the image) //
    ///////////////////////////////////////

    // Windows BITMAPINFOHEADER header

    // Read the bitmap information header  (DIB header)
    // 10 items * 4 bytes each
    bmp_header header;
    unsigned long dib_header_length = fread(&header, 1, 40, in);
    if (dib_header_length != 40) {
        return IO_ERR_FILE_TRUNC;
    }
    from_little_endian((int32_t *) &header, 10);

    // Check reported size of DIB header is 40 B (10 items * 4 bytes each)
    if (header.size != 40) {
        return IO_ERR_FILE_HEADER;
    }

    // Check no. of planes is 1
    if (header.planes != 1) {
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
        return IO_ERR_UNSUPPORTED;
    }

    // Check if compression type is supported by this library
    // Currently, we only work with uncompressed files.
    if (header.compression != 0) {
        return IO_ERR_UNSUPPORTED;
    }

    // TODO: Color palette

    puts("\nDIB header info");
    printf("Size of header (should be 40) = %d bytes\n", header.size);
    printf("Size of image (could be 0 if BI_RGB) = %d\n", header.image_size);
    printf("No. of components = %d\n", num_components);
    printf("Width and Height = %dx%d px\n", header.width, header.height);
    printf("Compression method = %d\n", header.compression);

    // Close file
    fclose(in);
    return 0;
}
