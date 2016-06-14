/* Copyright 2004, Boutell.Com, Inc. and Tobacco Documents Online.
   Released under the GNU General Public License, version 2 or later. */

#include <stdio.h>
#include <setjmp.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <tiffio.h>
#include <png.h>

void usage(char *s);
void die(char *s);

#define GetBWPixel(raster, y, x) \
	(raster[y * bwidth + (x >> 3)] & (1 << ((~x) & 7))) ? 0 : 255;

#define sign(a) ((a > 0) ? 1 : -1)

int main(int argc, char *argv[])
{
	int result;
	TIFF *tiff;
	u_long width, height;
	int bwidth;
	u_char *raster = 0;
	int r;
	int i;
	png_structp png_ptr;
	png_infop info_ptr;
	png_bytep row_pointer, row;
	int page = 0;
	int angle = 0;
	int nwidth = -1;
	int lr = 0;
	int tb = 0;
	int antialias = 0;
	int xisy = 0;
	int xisflipped = 0;
	int yisflipped = 0;
	uint16 bitsPerSample;
	char *tiffFile = 0;
	char *pngFile = 0;	
	FILE *out;
	int y;
	u_long num;
	u_long denom;
	char error[1024];
	for (i = 1; (i < argc); i++) {
		if (!strcmp(argv[i], "-p")) {
			if (argc == (i + 1)) {
				die("-p expects a page number");
			}
			page = atoi(argv[i + 1]) - 1;
			i++;
			if (page < 0) {
				die("-p expects a page number >= 1");
			}
		} else if (!strcmp(argv[i], "-r")) {
			if (argc == (i + 1)) {
				die("-r expects an angle");
			}
			angle = atoi(argv[i + 1]);
			i++;	
			if ((angle < 0) || (angle > 270) ||
				(angle % 90))
			{
				die("-a expects an angle of 0, 90, 180, or 270 degrees");
			}
		} else if (!strcmp(argv[i], "-w")) {
			if (argc == (i + 1)) {
				die("-w expects a width in pixels");
			}
			nwidth = atoi(argv[i + 1]);
			i++;
			if (nwidth <= 0) {
				die("-w expects a positive width in pixels");
			}
		} else if (!strcmp(argv[i], "-lr")) {
			lr = 1;
		} else if (!strcmp(argv[i], "-tb")) {		
			tb = 1;
		} else if (!strcmp(argv[i], "-a")) {		
			antialias = 1;
		} else if ((!tiffFile) && (argv[i][0] != '-')) {
			tiffFile = argv[i];
		} else if ((!pngFile) && (argv[i][0] != '-')) {
			pngFile = argv[i];
		} else {
			usage("unknown parameter");
		}
	}
	if (!tiffFile) {
		usage("tiff filename is required");
	}
	tiff = TIFFOpen(tiffFile, "rb");
	if (!tiff) {
		die("Can't open file");
	}	
	if (!TIFFSetDirectory(tiff, page)) {
		die("Can't access page number requested");
	}
	(void) TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
	if (nwidth == -1) {
		num = 1;
		denom = 1;
	} else {
		if (!xisy) {
			num = width;
			denom = nwidth;
		} else {
			num = height;
			denom = nwidth;
		}
	}	
	if ((!num) || (!denom)) {
		die("Width and height must both be nonzero");
	}
	(void) TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);
	(void) TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
	if (bitsPerSample != 1) {
		die("Sorry, only 1-bit (scan/fax) TIFFs are supported");
	}
	bwidth = (int) TIFFScanlineSize(tiff);
	raster = calloc(1, bwidth * height);
	if (!raster) {
		die("Memory allocation error");
	}
	for (y = 0; (y < height); y++) {
		TIFFReadScanline(tiff, raster + bwidth * y, y, 0);
	}
	switch (angle) {
		case 0:
		break;
		case 90:
		xisy = !xisy;
		yisflipped = !yisflipped;
		break;
		case 180:
		xisflipped = !xisflipped;
		yisflipped = !yisflipped;
		break;
		case 270:
		xisy = !xisy;
		xisflipped = !xisflipped;
		break;
	}
	if (lr) {
		xisflipped = !xisflipped;
	}
	if (tb) {
		yisflipped = !yisflipped;
	}
	row = calloc(sizeof(u_char), (xisy) ? height : width);
	if (!row) {
		die("Memory allocation error");
	}	
	png_ptr = png_create_write_struct(
		PNG_LIBPNG_VER_STRING, 0, 0, 0);		
	if (!png_ptr) {
		die("Cannot allocate png_structp");
	}
	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		die("Cannot allocate png_infop");
	}
	if (setjmp(png_jmpbuf(png_ptr))) {
		die("Error on png write");
	}
	if (pngFile) {
		out = fopen(pngFile, "wb");
		if (!out) {
			die("Cannot create output file");
		}
	} else {
		out = stdout;
		SET_BINARY(STDOUT_FILENO);		
	}
	png_init_io(png_ptr, out);
	/* Turning off filtering yields a large speed improvement at a 
		modest price in file size */
	png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
	if (nwidth == -1) {
		nwidth = width;
	}
	png_set_IHDR(png_ptr, info_ptr, ((xisy) ? height : width) * denom / num,
		((xisy) ? width : height) * denom / num,
		8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);
	for (y = 0; (y < info_ptr->height); y++) {
		int x;
		u_char *p;
		row_pointer = row;
		p = row;
		for (x = 0; (x < info_ptr->width); x++) {
			int tx, ty;
			int accum = 0, total = 0;
			int ty1 = (xisy ? x : y) * num / denom;
			int tx1 = (xisy ? y : x) * num / denom;
			int tx2, ty2;
			int xsteps, ysteps;
			if (!antialias) {
				tx2 = tx1 + 1;
				ty2 = ty1 + 1;
			} else {
				ty2 = (xisy ? (x + 1) : (y + 1)) * num / denom;
				tx2 = (xisy ? (y + 1) : (x + 1)) * num / denom;
			}
			ysteps = abs(ty2 - ty1);
			xsteps = abs(tx2 - tx1);
			if (xisflipped) {
				tx1 = width - 1 - tx1;
				tx2 = width - 1 - tx2;
			}	
			if (yisflipped) {
				ty1 = height - 1 - ty1;
				ty2 = height - 1 - ty2;
			}	
			ty = ty1;	
			while (ty != ty2) {
				tx = tx1;
				while (tx != tx2) {
					accum += GetBWPixel(raster, ty, tx);
					total++;
					tx += sign(tx2 - tx1);
				}
				ty += sign(ty2 - ty1);
			}
			if (total > 0) {
				*p = accum / total;	
			}
			p++;
		}
		png_write_row(png_ptr, row_pointer);	
		
	}
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	if (out != stdout) {
		fclose(out);
	}
	return 0;
}

void usage(char *s)
{
	fprintf(stderr, "%s\n\n", s);	
	die("fax2png, version 1.0, copyright 2004, Boutell.Com, Inc. and\n"
		"Tobacco Documents Online.\n\n"
		"Usage: fax2png tifffilename [pngfilename] [-p pagenumber] [-w width]\n"
		"	[-r 0|90|180|270] [-lr] [-tb] [-a] (the -a option is strongly\n"
		"	recommended for attractive results)");
}

void die(char *s)
{
	fprintf(stderr, "%s\n", s);
	exit(1);
}

