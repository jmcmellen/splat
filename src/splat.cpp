/****************************************************************************\
 *	   SPLAT!: An RF Signal Path Loss And Terrain Analysis Tool		  *
 ******************************************************************************
 *		 Project started in 1997 by John A. Magliacane, KD2BD 		 *
 *			  Last update: 31-May-2019			                             *
 ******************************************************************************
 *		 Please consult the documentation for a complete list of		 *
 *		 individuals who have contributed to this project. 			 *
 ******************************************************************************
 *																			*
 *  This program is free software; you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the	 *
 *  Free Software Foundation; either version 2 of the License or any later	*
 *  version.									 *
 * 										 *
 *  This program is distributed in the hope that it will useful, but WITHOUT  *
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or	 *
 *  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License	 *
 *  for more details.								 *
 *										 *
 \****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#ifndef _WIN32
#include <unistd.h>
#include <limits.h>
#include <libgen.h> /* dirname() and basename() */
#else
#define unlink _unlink
#endif

#include "bzlib.h"
#include "zip.h"

#define HAVE_LIBPNG
#ifndef _WIN32
#define HAVE_LIBJPEG
#endif

#ifdef HAVE_LIBPNG
#include "png.h"
#endif
#ifdef HAVE_LIBJPEG
#include "jpeglib.h"
#endif

#include "fontdata.h"
#include "workqueue.hpp"

#ifdef _WIN32
#define PATHSEP '\\'
#define strdup _strdup
#else
#define PATHSEP '/'
#endif

#define GAMMA 2.5
#define BZBUFFER 65536

#define SPLAT_VERSION "2.0.0b4"
#define SPLAT_NAME "SPLAT! MT"

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef TWOPI
#define TWOPI 6.283185307179586
#endif

#ifndef HALFPI
#define HALFPI 1.570796326794896
#endif

#define DEG2RAD 1.74532925199e-02
#define EARTHRADIUS 20902230.97
#define	METERS_PER_MILE 1609.344
#define METERS_PER_FOOT 0.3048
#define	KM_PER_MILE 1.609344
#define FOUR_THIRDS 1.3333333333333

#define MAX_PATH_LEN 255
#define MAX_LINE_LEN 128
#define MAX_SITE_NAME_LEN 64

#ifndef min
#define min(i, j) ( i < j ? i : j)
#endif
#ifndef max
#define max(i, j) ( i > j ? i : j)
#endif

#define DO_KMZ 1

/* avg must be a float or double or you get rounding down because of truncation
 * ignores values if n == 0 to avoid a divide-by-zero error
 */
#define UPDATE_RUNNING_AVG(avg, latest, n) do { \
    (n) && ((avg) = (avg) + ((float)((latest) - (avg)))/(float)(n)); \
} while (0)

static double avgpathlen = 0.0;
static int totalpaths = 0;

/*****************************
 * Typedefs
 *****************************/

typedef struct Site {
    double lat;
    double lon;
    double alt;
    unsigned char amsl_flag;
    char name[MAX_SITE_NAME_LEN];
    char filename[MAX_PATH_LEN];
} Site;

/* data about a path between two points
 * Note, Mac OSX has a default pthread stack size of 512KB, so trying to
 * allocate one of these on the stack will cause it to blow up with a
 * permissions (access) error. Allocating it on the heap via malloc() or
 * new() is, of course, OK.
 */
typedef struct Path {
    unsigned long arysize;
    double *lat;
    double *lon;
    double *elevation;
    double *distance;
    int length;
} Path;

// digital elevation model data
typedef struct DEM {
    int min_north;
    int max_north;
    int min_west;
    int max_west;
    int max_el;
    int min_el;
    unsigned long arysize;
    short **data;
    unsigned char **mask;
    unsigned char **signal;
} DEM;

#define ANTENNA_RADIALS 361
#define ANTENNA_ANGLES 1001
#define NO_ANTENNA_DATA (-1.0)
typedef struct LongleyRiceData {
    double eps_dielect;
    double sgm_conductivity;
    double eno_ns_surfref;
    double frq_mhz;
    double conf;
    double rel;
    double erp;
    int radio_climate;
    int pol;
    double antenna_pattern[ANTENNA_RADIALS][ANTENNA_ANGLES];
} LongleyRice;

typedef struct Region {
    unsigned char color[32][3];
    int level[32];
    int levels;
} Region;

typedef enum ImageType {
    IMAGETYPE_PPM = 0,
#ifdef HAVE_LIBPNG
    IMAGETYPE_PNG,
#endif
#ifdef HAVE_LIBJPEG
    IMAGETYPE_JPG,
#endif
} ImageType;

typedef enum SDFCompressType {
    SDF_COMPRESSTYPE_NONE = 0,
    SDF_COMPRESSTYPE_BZIP2
} SDFCompressType;

struct SDFCompressFormat {
    SDFCompressType type;
    const char* suffix;
};


/*****************************
 * Globals
 *****************************/
typedef enum AppMode {
    APPMODE_NORMAL = 0,
    APPMODE_HD = 1
} AppMode;
AppMode appmode = APPMODE_NORMAL;

#define MAXPAGESIDES 8
const int appArraySize[2][MAXPAGESIDES] = {
    { 4950, 4950, 10870, 19240, 30025, 43217, 58813, 76810 },  /* Normal mode */
    { 5092, 14844, 32600, 57713, 90072, 129650, 176437, 230430 } /* HD mode */
};

char sdf_path[MAX_PATH_LEN], home_sdf_path[MAX_PATH_LEN], dashes[80];

double	earthradius, max_range=0.0, forced_erp=-1.0, dpp, ppd,
        fzone_clearance=0.6, forced_freq, clutter;

int	min_north=90, max_north=-90, min_west=360, max_west=-1;
int max_elevation=-32768, min_elevation=32768;
int ippd, mpi, maxpagesides = MAXPAGESIDES;
int bzerror, contour_threshold;

bool got_elevation_pattern = false, got_azimuth_pattern = false;
bool metric = false, dbm = false, smooth_contours = false;
bool gpsav=false, itwom=false;

int verbose=1;

bool opened=false; /* BZip file flag */

LongleyRiceData LR;
Region region;

DEM **aDEM = NULL;
int demCount = 0;

/****************************
 * Color handling stuff
 ****************************/

/* 100: area plot image size is about 550KB (best quality)
 *  95: area plot image size is about 290KB.
 *  90: area plot image size is about 200KB. This is about what "convert"'s default is.
 *  80: area plot image size is about 130KB, but has noticeable artifacts when zoomed in.
 */
#define DEFAULT_JPEG_QUALITY 90

typedef uint32_t Pixel;

#ifndef _WIN32
#define RGB(r,g,b) ( ((uint32_t)(uint8_t)r)|((uint32_t)((uint8_t)g)<<8)|((uint32_t)((uint8_t)b)<<16) )

#define GetRValue(RGBColor) (uint8_t) (RGBColor & 0xff)
#define GetGValue(RGBColor) (uint8_t) ((((uint32_t)(RGBColor)) >> 8) & 0xff)
#define GetBValue(RGBColor) (uint8_t) ((((uint32_t)(RGBColor)) >> 16) & 0xff)
#endif

#define COLOR_RED				RGB(255, 0, 0)
#define COLOR_LIGHTCYAN			RGB(128,128,255)
#define COLOR_GREEN				RGB(0,255,0)
#define COLOR_DARKGREEN			RGB(0,100,0)
#define COLOR_DARKSEAGREEN1		RGB(193,255,193)
#define COLOR_CYAN				RGB(0,255,255)
#define COLOR_YELLOW			RGB(255,255,0)
#define COLOR_GREENYELLOW		RGB(173,255,47)
#define COLOR_MEDIUMSPRINGGREEN	RGB(0,250,154)
#define COLOR_MEDIUMVIOLET		RGB(147,112,219)
#define COLOR_PINK				RGB(255,192,203)
#define COLOR_ORANGE			RGB(255,165,0)
#define COLOR_SIENNA			RGB(255,130,71)
#define COLOR_BLANCHEDALMOND	RGB(255,235,205)
#define COLOR_DARKTURQUOISE		RGB(0,206,209)
#define COLOR_TAN				RGB(210,180,140)
#define COLOR_GOLD2				RGB(238,201,0)
#define COLOR_MEDIUMBLUE		RGB(0,0,170)
#define COLOR_WHITE				RGB(255,255,255)
#define COLOR_BLACK				RGB(0,0,0)




/*****************************
 * External stuff for ITWOM. I don't know why we don't just include a .h file.
 *****************************/
#ifdef ITM_ELEV_DOUBLE
#define elev_t double
#else
#define elev_t float
#endif

#ifdef __cplusplus
extern "C" {
#endif
    void point_to_point(elev_t elev[], double tht_m, double rht_m,
                        double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
                        double frq_mhz, int radio_climate, int pol, double conf,
                        double rel, double *dbloss, char *strmode, int *errnum);

    void point_to_point_ITM(elev_t elev[], double tht_m, double rht_m,
                        double eps_dielect, double sgm_conductivity, double eno_ns_surfref,
                        double frq_mhz, int radio_climate, int pol, double conf,
                        double rel, double *dbloss, char *strmode, int *errnum);

    double ITWOMVersion();
#ifdef __cplusplus
}
#endif

/*****************************
 * System-level utilities
 *****************************/

#ifndef _WIN32
#include <unistd.h>

unsigned long long getTotalSystemMemory()
{
    unsigned long pages = sysconf(_SC_PHYS_PAGES);
    unsigned long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
#else
#include <windows.h>

unsigned long long getTotalSystemMemory()
{
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
}

#endif

/*****************************
 * Image writing functions. This should just be a class.
 *****************************/

typedef struct ImageWriter_st
{
    bool initialized;

    FILE *fp;
    ImageType imagetype;
    int width;
    int height;

    int xoffset;

    unsigned char *imgline;

#ifdef HAVE_LIBPNG
    png_structp png_ptr;
    png_infop info_ptr;
#endif
#ifdef HAVE_LIBJPEG
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
#endif
} ImageWriter;

int ImageWriterInit(ImageWriter *iw, const char* filename, ImageType imagetype, int width, int height)
{
    memset(iw, 0, sizeof(ImageWriter));

    if ( (iw->fp=fopen(filename,"wb")) == NULL) {
        return -1;
    }

    iw->imagetype = imagetype;
    iw->width = width;
    iw->height = height;
    iw->xoffset = 0;

    iw->imgline = (unsigned char*)malloc(3 * width);

    // XXX TODO: error handling
    switch (iw->imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            iw->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
            iw->info_ptr = png_create_info_struct(iw->png_ptr);
            png_init_io(iw->png_ptr, iw->fp);
            png_set_IHDR(iw->png_ptr, iw->info_ptr,
                         iw->width, iw->height,
                         8, /* 8 bits per color or 24 bits per pixel for RGB */
                         PNG_COLOR_TYPE_RGB,
                         PNG_INTERLACE_NONE,
                         PNG_COMPRESSION_TYPE_DEFAULT,
                         PNG_FILTER_TYPE_BASE);
            png_set_compression_level(iw->png_ptr, 6);  /* default is Z_DEFAULT_COMPRESSION; see zlib.h */
            png_write_info(iw->png_ptr, iw->info_ptr);
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            iw->cinfo.err = jpeg_std_error(&iw->jerr);
            jpeg_create_compress(&iw->cinfo);
            jpeg_stdio_dest(&iw->cinfo, iw->fp);
            iw->cinfo.image_width = iw->width;
            iw->cinfo.image_height = iw->height;
            iw->cinfo.input_components = 3;	   /* # of color components per pixel */
            iw->cinfo.in_color_space = JCS_RGB;   /* colorspace of input image */
            jpeg_set_defaults(&iw->cinfo); /* default compression */
            jpeg_set_quality(&iw->cinfo, DEFAULT_JPEG_QUALITY, TRUE); /* possible range is 0-100 */
            jpeg_start_compress(&iw->cinfo, TRUE); /* start compressor. */
            break;
#endif
        case IMAGETYPE_PPM:
            fprintf(iw->fp,"P6\n%u %u\n255\n",iw->width,iw->height);
    }

    iw->initialized = true;
    return 0;
}

void ImageWriterAppendPixel(ImageWriter *iw, Pixel pixel)
{
    if (!iw->initialized)
        return;

    if ((iw->xoffset+3) > (iw->width*3)) {
        return;
    }

    iw->imgline[iw->xoffset++]=GetRValue(pixel);
    iw->imgline[iw->xoffset++]=GetGValue(pixel);
    iw->imgline[iw->xoffset++]=GetBValue(pixel);
}

void ImageWriterEmitLine(ImageWriter *iw)
{
    if (!iw->initialized)
        return;

    switch (iw->imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            png_write_row(iw->png_ptr, (png_bytep)(iw->imgline));
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            jpeg_write_scanlines(&iw->cinfo, &iw->imgline, 1);
            break;
#endif
        case IMAGETYPE_PPM:
            fwrite(iw->imgline, 3, iw->width, iw->fp);
            break;
    }
    iw->xoffset = 0;
}

void ImageWriterFinish(ImageWriter *iw)
{
    if (!iw->initialized)
        return;

    switch (iw->imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            png_write_end(iw->png_ptr, iw->info_ptr);
            png_destroy_write_struct(&iw->png_ptr, &iw->info_ptr);
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            jpeg_finish_compress(&iw->cinfo);
            jpeg_destroy_compress(&iw->cinfo);
            break;
#endif
        case IMAGETYPE_PPM:
            /* do nothing */
            ;
    }

    fclose(iw->fp);
    free(iw->imgline);
}

/*****************************
 * Functions for allocating Paths
 ****************************/
int AppArraySize(AppMode mode)
{
    if (mode == APPMODE_HD) {
        return 3600;
    }
    return 1200;
}

Path *AllocatePath()
{
    Path *path = (Path*)malloc(sizeof(Path));
    if (!path) {
        return NULL;
    }
    memset(path, 0, sizeof(Path));
    path->arysize = appArraySize[appmode][maxpagesides-1];

    path->lat=(double*)malloc(path->arysize*sizeof(double));
    path->lon=(double*)malloc(path->arysize*sizeof(double));
    path->elevation=(double*)malloc(path->arysize*sizeof(double));
    path->distance=(double*)malloc(path->arysize*sizeof(double));

    if (path->lat && path->lon && path->elevation && path->distance) {
        return path;
    }

    /* if we got here, there was an error */
    if (path->lat) free(path->lat);
    if (path->lon) free(path->lon);
    if (path->elevation) free(path->elevation);
    if (path->distance) free(path->distance);

    free(path);

    return NULL;
}

void DestroyPath(Path *path)
{
    if (path->lat) free(path->lat);
    if (path->lon) free(path->lon);
    if (path->elevation) free(path->elevation);
    if (path->distance) free(path->distance);
    free(path);
}

unsigned long SizeofPath(AppMode mode, int pagesides)
{
    unsigned long arysize = appArraySize[mode][pagesides-1];
    return sizeof(Path) + (sizeof(double)*arysize*4);
}

/*****************************
 * Utility functions
 *****************************/
char* basename_s(char* path) {
    if (!path) return (char*)""; /* const string */

    int i = (int)strlen(path);
    for ( ; i >= 0 && i != '/' && i != '\\'; --i); /* walk backwards until we find a path separator */
    ++i; /* i is either -1 or the index of the slash, so move forward */
    return &(path[i]);
}

/* Strip the extension off a path
 *    path: path to strip
 *     ext: buffer to hold the extension
 *  extlen: size of ext INCLUDING the null byte
 */
void stripExtension(char *path, char* ext, size_t extlen)
{
    char *p;
    size_t len;

    p = strrchr(path, '.');

    if (p == NULL) {
        ext[0] = '\0';
        return;
    }

    *p++ = '\0'; /* chop off extension and move ahead */

    len = min(strlen(p)+1, extlen);
    memcpy(ext, p, len);
    ext[len] = '\0';
}

#ifdef _WIN32
/* gnuplot's internal interpreter doesn't handle single backslashes
 * properly on windows. Our options are to either escape them by turning
 * them into double backslashes or converting them to forward slashes.
 * Here we do an in-place conversion to forward slashes so that we don't
 * have to reallocate any memory. */
void convertBackslashes(char *path)
{
    size_t i;
    for (i=strlen(path)-1; i; --i) {
        if (path[i] == '\\') path[i]='/';
    }
}
#else
/* empty macro */
#define convertBackslashes(x) 
#endif

#ifndef _WIN32
/* Windows has a useful function called _splitpath_s(). This is
 * intended to replicate that. */
int _splitpath_s(
                 const char * path,
                 char * dir,
                 size_t dirNumberOfElements,
                 char * fname,
                 size_t nameNumberOfElements,
                 char * ext,
                 size_t extNumberOfElements
                )
{
    char *p;

    if (path == NULL) {
        return EINVAL;
    }

    char *fullpath = strdup(path);

    /* Peel off everything up to (but not including) the final slash.
     * If there is no slash (local file), returns an empty string.
     * Guaranteed to be NULL terminated.
     */
    if (dir != NULL) {
        char *tmp = strdup(fullpath); /* dirname munges its arg so copy*/
        char *dname = dirname(tmp); /* If there is no slash, gets ".". */
        if (dirNumberOfElements < strlen(dname)+1) {
            free(tmp);
            free(fullpath);
            return ERANGE;
        }
        if (dname[0] == '.' && dname[1] == '\0') {
            dir[0] = '\0';
        }
        else {
            memcpy(dir, dname, strlen(dname) + 1);
        }
        free(tmp);
    }

    /* Everything after the last slash, or the whole thing if
     * there is no slash. guaranteed to be null-terminated.
     */
    char *base = basename(fullpath);

    if ( ((p = strchr(base, '.')) == NULL) || (strlen(p) < 2)) {
        /* Easy case: no suffix */
        if (ext) {
            if (extNumberOfElements < 1) {
                free(fullpath);
                return ERANGE;
            }
            ext[0] = '\0';
        }
    } else {
        /* peel off everything from the first dot to the end. Note
         * that multiple suffixes get treated as one; i.e. "foo.tar.gz"
         * gets split into "foo" and "tar.gz".
         */
        *p++ = '\0';
        if (ext) {
            if (extNumberOfElements < strlen(p)) {
                free(fullpath);
                return ERANGE;
            }
            memcpy(ext, p, strlen(p) + 1);
        }
    }

    if (fname) {
        if (nameNumberOfElements < strlen(base)+1) {
            free(fullpath);
            return ERANGE;
        }
        memcpy(fname, base, strlen(base) + 1);
    }

    free(fullpath);
    return 0;
}
#endif

/* Copy a filename, up to but not including the suffix.
 * If newSuffix is specified, add that on.
 * The new filename is malloced and must be freed when
 * the caller is done with it.
 */
char* copyFilename(char *fname, const char *newSuffix) {
    char pathbuf[1027], filebuf[512], suffixbuf[64];
    char *newFname;

#ifdef _WIN32
    char drivebuf[3];
    drivebuf[0] = '\0';
    if (_splitpath_s(fname, drivebuf, 3, pathbuf, 1024, filebuf, 512, suffixbuf, 64) != 0) {
        return NULL;
    }
    size_t offset = strlen(drivebuf); /* includes room for colon */
    if ((offset > 0) && (strlen(pathbuf) > 0)) {
        memmove(pathbuf+offset, pathbuf, strlen(pathbuf));
        memcpy(pathbuf, drivebuf, offset);
    }
#else
    if (_splitpath_s(fname, pathbuf, 1024, filebuf, 512, suffixbuf, 64) != 0) {
        return NULL;
    }
#endif

    if (newSuffix && (strlen(newSuffix) > 0)) {
        size_t len = strlen(pathbuf) + strlen(filebuf) + strlen(newSuffix) + 3;
        newFname = (char*)malloc(len);
        if (strlen(pathbuf)==0 || strcmp(pathbuf, ".")==0) {  /* if it's "./foo", just make it "foo" */
            snprintf(newFname, len, "%s.%s", filebuf, newSuffix);
        } else {
            snprintf(newFname, len, "%s%c%s.%s", pathbuf, PATHSEP, filebuf, newSuffix);
        }
    } else {
        size_t len = strlen(pathbuf) + strlen(filebuf) + 2;
        newFname = (char*)malloc(len);
        if (strlen(pathbuf) == 0 || strcmp(pathbuf, ".")==0) {  /* if it's "./foo", just make it "foo" */
            snprintf(newFname, len, "%s", filebuf);
        } else {
            snprintf(newFname, len, "%s%c%s", pathbuf, PATHSEP, filebuf);
        }
    }
    return newFname;
}


/* Perform linear interpolation between quantized contour
 * levels displayed in field strength and path loss maps.
 * If signal level x0 corresponds to color level y0, signal
 * level x1 corresponds to color level y1, and signal level
 * n falls somewhere between x0 and x1, determine what
 * color value n corresponds to between y0 and y1.
 */
int interpolate(int y0, int y1, int x0, int x1, int n)
{
    int result=0;
    double delta_x, delta_y;

    if (n<=x0)
        return y0;

    if (n>=x1)
        return y1;

    if (y0==y1)
        return y0;

    if (x0==x1)
        return y0;

    delta_y=(double)(y0-y1);
    delta_x=(double)(x0-x1);

    result=y0+(int)ceil((delta_y/delta_x)*(n-x0));

    return result;
}

/* Implements the arc cosine function, returning a value
 * between 0 and TWOPI.
 */
double arccos(double x, double y)
{
    double result=0.0;

    if (y>0.0)
        result=acos(x/y);

    if (y<0.0)
        result=PI+acos(x/y);

    return result;
}

/* Normalizes the argument to an integer angle
 * between 0 and 180 degrees.
 */
int ReduceAngle(double angle)
{
    double temp;

    temp=acos(cos(angle*DEG2RAD));

    return (int)rint(temp/DEG2RAD);
}

/* This function returns the short path longitudinal
 * difference between longitude1 and longitude2 
 * as an angle between -180.0 and +180.0 degrees.
 * If lon1 is west of lon2, the result is positive.
 * If lon1 is east of lon2, the result is negative.
 */
double LonDiff(double lon1, double lon2)
{
    double diff;

    diff=lon1-lon2;

    if (diff<=-180.0)
        diff+=360.0;

    if (diff>=180.0)
        diff-=360.0;

    return diff;
}

/* Converts decimal degrees to degrees, minutes, seconds,
 * (DMS) and returns the result as a character string.
 * Uses an internal static buffer; NOT THREADSAFE.
 */
char *dec2dms(double decimal)
{
    char	sign;
    int	degrees, minutes, seconds;
    double	a, b, c, d;

    static char buf[MAX_LINE_LEN];

    if (decimal<0.0)
    {
        decimal=-decimal;
        sign=-1;
    }

    else
        sign=1;

    a=floor(decimal);
    b=60.0*(decimal-a);
    c=floor(b);
    d=60.0*(b-c);

    degrees=(int)a;
    minutes=(int)c;
    seconds=(int)d;

    if (seconds<0)
        seconds=0;

    if (seconds>59)
        seconds=59;

    buf[0]=0;
    snprintf(buf,MAX_LINE_LEN,"%d%c %d\' %d\"", degrees*sign, 176, minutes, seconds);
    return (buf);
}

/*****************************
 * Functions for allocating and manipulating DEMs
 *****************************/

/* Allocate a DEM struct on the heap.
 */
DEM *AllocateDEM()
{
    /* not allowed to allocate any more */
    if (demCount >= (maxpagesides*maxpagesides)) {
        return NULL;
    }

    DEM *dem = (DEM*)malloc(sizeof(DEM));
    if (!dem) {
        return NULL;
    }
    memset(dem, 0, sizeof(DEM));
    dem->arysize = ippd;
    if (dem->arysize < 1) dem->arysize = 1;

    do {
        dem->data=(short**)malloc(dem->arysize*sizeof(short*));
        dem->mask=(unsigned char**)malloc(dem->arysize*sizeof(unsigned char*));
        dem->signal=(unsigned char**)malloc(dem->arysize*sizeof(unsigned char*));
        if (!dem->data || !dem->mask || !dem->signal) {
            break;
        }

        for (unsigned long i = 0; i<dem->arysize; ++i) {
            dem->data[i]=(short*)malloc(dem->arysize*sizeof(short));
            dem->mask[i]=(unsigned char*)malloc(dem->arysize*sizeof(unsigned char));
            dem->signal[i]=(unsigned char*)malloc(dem->arysize*sizeof(unsigned char));
            // XXX TODO: error check these allocations
        }

        dem->min_el=32768;
        dem->max_el=-32768;
        dem->min_north=90;
        dem->max_north=-90;
        dem->min_west=360;
        dem->max_west=-1;

        return dem;

    } while (0);

    /* if we got here, there was an error */
    if (dem->signal) free(dem->signal);
    if (dem->mask) free(dem->mask);
    if (dem->data) free(dem->data);
    free(dem);

    return NULL;
}

void DestroyDEM(DEM *dem)
{
    for (unsigned long i = 0; i<dem->arysize; ++i) {
        free(dem->data[i]);
        free(dem->mask[i]);
        free(dem->signal[i]);
    }
    if (dem->signal) free(dem->signal);
    if (dem->mask) free(dem->mask);
    if (dem->data) free(dem->data);
    free(dem);
}

unsigned long SizeofDEM(AppMode mode)
{
    unsigned long arysize = AppArraySize(mode);
    return sizeof(DEM) + 
        (sizeof(short*)*arysize) +
        (sizeof(short)*arysize*arysize) + 
        (sizeof(unsigned char*)*arysize*2) +
        (sizeof(unsigned char)*arysize*arysize*2);
}

int InitDEMs()
{
    aDEM = (DEM**)calloc(maxpagesides*maxpagesides, sizeof(DEM*));
    if (aDEM) { return 0; };
    return -1;  // Insufficient memory?
}

void FreeDEMs()
{
    for (int i=0; i< demCount; ++i) {
        DestroyDEM(aDEM[i]);
    }
    free(aDEM);
    aDEM = NULL;
}

/* Appends a DEM to the global adem[] array.u
 */
int AppendDEM(DEM *dem)
{
    if (aDEM == NULL || (demCount > (maxpagesides*maxpagesides)) ) {
        return -1;
    }
    aDEM[demCount++] = dem;
    return 0;
}

/* Returns a DEM whose boundaries are exactly on the given boundaries,
 * or NULL if not found.
 */
DEM *FindDEM_Explicit(double minlat, double minlon, double maxlat, double maxlon)
{
    for (int i=0; i<demCount; ++i) {
        if (minlat==aDEM[i]->min_north &&
            minlon==aDEM[i]->min_west &&
            maxlat==aDEM[i]->max_north &&
            maxlon==aDEM[i]->max_west)
            return aDEM[i];
    }
    return NULL;
}

/* Returns the first DEM containing the lat/long,
 * or -1 if not found.
 *
 * x and y will contain the offsets into the DEM array of
 * the coordinate.
 */
DEM *FindDEM(double lat, double lon, int &x, int &y)
{
    for (int i=0; i<demCount; ++i)
    {
        x=(int)rint(ppd*(lat-aDEM[i]->min_north));
        y=mpi-(int)rint(ppd*(LonDiff(aDEM[i]->max_west,lon)));

        if (x>=0 && x<=mpi && y>=0 && y<=mpi)
            return aDEM[i];
    }
    return NULL;
}


/* Lines, text, markings, and coverage areas are stored in a
 * mask that is combined with topology data when topographic
 * maps are generated by SPLAT!.  This function sets and resets
 * bits in the mask based on the latitude and longitude of the
 * area pointed to.
 */
int PutMask(double lat, double lon, int value)
{
    int	x, y;
    DEM *dem; 

    dem = FindDEM(lat, lon, x, y);
    if (!dem)
        return -1;

    dem->mask[x][y]=value;
    return ((int)dem->mask[x][y]);
}

/* Lines, text, markings, and coverage areas are stored in a
 * mask that is combined with topology data when topographic
 * maps are generated by SPLAT!.  This function sets bits in
 * the mask based on the latitude and longitude of the area
 * pointed to.
 */
int OrMask(double lat, double lon, int value)
{
    int	x, y;
    DEM *dem;

    dem = FindDEM(lat, lon, x, y);
    if (!dem)
        return -1;

    dem->mask[x][y]|=value;
    return ((int)dem->mask[x][y]);
}

/* Returns the mask bits based on the latitude and
 * longitude given.
 */
int GetMask(double lat, double lon)
{
    return (OrMask(lat,lon,0));
}

/* Writes a signal level (0-255) at the specified location
 * lat/long location in the DEM for later recall.
 * Returns the signal value just set.
 */
int PutSignal(double lat, double lon, unsigned char signal)
{
    int	x, y;
    DEM *dem;

    dem = FindDEM(lat, lon, x, y);
    if (!dem)
        return 0;

    dem->signal[x][y]=signal;
    return (dem->signal[x][y]);
}

/* This function reads the signal level (0-255) at the
 *  specified location that was previously written by the
 * complimentary PutSignal() function.
 */
unsigned char GetSignal(double lat, double lon)
{
    int	x, y;
    DEM *dem;

    dem = FindDEM(lat, lon, x, y);
    if (!dem)
        return 0;

    return (dem->signal[x][y]);
}

/* This function returns the elevation (in feet) of any location
 * represented by the digital elevation model data in memory.
 * Function returns -5000.0 for locations not found in memory.
 */
double GetElevation(Site location)
{
    int	x, y;
    DEM *dem;

    dem = FindDEM(location.lat, location.lon, x, y);
    if (!dem)
        return -5000.0;

    return (3.28084*(double)(dem->data[x][y]));
}

/* This function adds a user-defined terrain feature
 * (in meters AGL) to the digital elevation model data
 * in memory.  Does nothing and returns 0 for locations
 * not found in memory.
 */
int AddElevation(double lat, double lon, double height)
{
    int	x, y;
    DEM *dem;

    dem = FindDEM(lat, lon, x, y);
    if (!dem)
        return 0;

    dem->data[x][y]+=(short)rint(height);
    return 1;
}

/* Sets the terrain elevation to the specified height above
 * mean sea level. Useful for compensating for canopy height
 * (ie: buildings and trees) included in SRTM and ASTER
 * topographic data sets by setting the terrain height to zero
 * when the height of an antenna above mean sea level is known.
 * Returns 0 for locations not found in memory.
 */
int SetElevation(double lat, double lon, double height)
{
    int x, y;
    DEM *dem;

    dem = FindDEM(lat, lon, x, y);
    if (!dem)
        return 0;

    dem->data[x][y]=(short)rint(height);
    return 1;
}


/*****************************
 * More utility functions, this time for Path
 *****************************/

/* Returns the great circle distance in miles between any
 * two site locations.
 */
double Distance(Site site1, Site site2)
{
    double	lat1, lon1, lat2, lon2, angle, distance;

    lat1=site1.lat*DEG2RAD;
    lon1=site1.lon*DEG2RAD;
    lat2=site2.lat*DEG2RAD;
    lon2=site2.lon*DEG2RAD;

    angle = sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos((lon1)-(lon2));
    angle = max(-1.0, min(1.0, angle));
    distance=3959.0*acos(angle);

    return distance;
}

/* This function returns the azimuth (in degrees) to the
 * destination as seen from the location of the source.
 */
double Azimuth(Site source, Site destination)
{
    double	dest_lat, dest_lon, src_lat, src_lon,
            beta, azimuth, diff, num, den, fraction;
    double angle;

    dest_lat=destination.lat*DEG2RAD;
    dest_lon=destination.lon*DEG2RAD;

    src_lat=source.lat*DEG2RAD;
    src_lon=source.lon*DEG2RAD;

    /* Calculate Surface Distance */

    angle = sin(src_lat)*sin(dest_lat)+cos(src_lat)*cos(dest_lat)*cos(src_lon-dest_lon);
    angle = max(-1.0, min(1.0, angle));
    beta=acos(angle);

    /* Calculate Azimuth */

    num=sin(dest_lat)-(sin(src_lat)*cos(beta));
    den=cos(src_lat)*sin(beta);
    fraction=num/den;

    /* Trap potential problems in acos() due to rounding */

    if (fraction>=1.0)
        fraction=1.0;

    if (fraction<=-1.0)
        fraction=-1.0;

    /* Calculate azimuth */

    azimuth=acos(fraction);

    /* Reference it to True North */

    diff=dest_lon-src_lon;

    if (diff<=-PI)
        diff+=TWOPI;

    if (diff>=PI)
        diff-=TWOPI;

    if (diff>0.0)
        azimuth=TWOPI-azimuth;

    return (azimuth/DEG2RAD);
}

/* Returns the angle of elevation (in degrees) of the destination
 * as seen from the source location. A positive result represents
 * an angle of elevation (uptilt), while a negative result represents
 * an angle of depression (downtilt), as referenced to a normal to
 * the center of the earth.
 */
double ElevationAngle(Site source, Site destination)
{
    double a, b, dx;

    a=GetElevation(destination)+destination.alt+earthradius;
    b=GetElevation(source)+source.alt+earthradius;

    dx=5280.0*Distance(source,destination);

    /* Apply the Law of Cosines */

    return ((180.0*(acos(((b*b)+(dx*dx)-(a*a))/(2.0*b*dx)))/PI)-90.0);
}

/*****************************
 * The one-and-only Path reading function
 *****************************/

/* Generates a sequence of latitude and longitude positions between
 * source and destination locations along a great circle path, and
 * stores elevation and distance information for points along that
 * path in the "path" structure.
 */
void ReadPath(Site &source, Site &destination, Path *path)
{
    int	c;
    double	azimuth, distance, lat1, lon1, beta, den, num,
            lat2, lon2, total_distance, dx, dy, path_length,
            miles_per_sample, samples_per_radian=68755.0;
    Site tempsite;

    lat1=source.lat*DEG2RAD;
    lon1=source.lon*DEG2RAD;

    lat2=destination.lat*DEG2RAD;
    lon2=destination.lon*DEG2RAD;

    if (ppd==1200.0)
        samples_per_radian=68755.0;

    if (ppd==3600.0)
        samples_per_radian=206265.0;

    azimuth=Azimuth(source,destination)*DEG2RAD;

    total_distance=Distance(source,destination);

    if (total_distance>(30.0/ppd))		/* > 0.5 pixel distance */
    {
        dx=samples_per_radian*acos(cos(lon1-lon2));
        dy=samples_per_radian*acos(cos(lat1-lat2));

        path_length=sqrt((dx*dx)+(dy*dy));		/* Total number of samples */

        miles_per_sample=total_distance/path_length;	/* Miles per sample */
    }

    else
    {
        c=0;
        dx=0.0;
        dy=0.0;
        path_length=0.0;
        miles_per_sample=0.0;
        total_distance=0.0;

        lat1=lat1/DEG2RAD;
        lon1=lon1/DEG2RAD;

        path->lat[c]=lat1;
        path->lon[c]=lon1;
        path->elevation[c]=GetElevation(source);
        path->distance[c]=0.0;
    }

    for (distance=0.0, c=0; (total_distance!=0.0 && distance<=total_distance && c<(int)path->arysize); c++, distance=miles_per_sample*(double)c)
    {
        beta=distance/3959.0;
        lat2=asin(sin(lat1)*cos(beta)+cos(azimuth)*sin(beta)*cos(lat1));
        num=cos(beta)-(sin(lat1)*sin(lat2));
        den=cos(lat1)*cos(lat2);

        if (azimuth==0.0 && (beta>HALFPI-lat1))
            lon2=lon1+PI;

        else if (azimuth==HALFPI && (beta>HALFPI+lat1))
            lon2=lon1+PI;

        else if (fabs(num/den)>1.0)
            lon2=lon1;

        else
        {
            if ((PI-azimuth)>=0.0)
                lon2=lon1-arccos(num,den);
            else
                lon2=lon1+arccos(num,den);
        }

        if (lon2<0.0)
        {
            while (lon2<0.0)
                lon2+=TWOPI;
        }

        if (lon2>TWOPI)
        {
            while (lon2>TWOPI)
                lon2-=TWOPI;
        }

        lat2=lat2/DEG2RAD;
        lon2=lon2/DEG2RAD;

        path->lat[c]=lat2;
        path->lon[c]=lon2;
        tempsite.lat=lat2;
        tempsite.lon=lon2;

        path->elevation[c]=GetElevation(tempsite);
        path->distance[c]=distance;
    }

    /* Make sure exact destination point is recorded at path.length-1 */

    if (c<(int)path->arysize)
    {
        path->lat[c]=destination.lat;
        path->lon[c]=destination.lon;
        path->elevation[c]=GetElevation(destination);
        path->distance[c]=total_distance;
        c++;
    }

    if (c<(int)path->arysize)
        path->length=c;
    else
        path->length=path->arysize-1;
}

/* Returns the angle of elevation (in degrees) of the destination as
 * seen from the source location, UNLESS the path between the sites is
 * obstructed, in which case it returns the elevation angle to the first
 * obstruction.
 * "er" represents the earth radius.
 */
double ElevationAngle2(Site source, Site destination, double er)
{
    int	x;
    char	block=0;
    double	source_alt, destination_alt, cos_xmtr_angle,
            cos_test_angle, test_alt, elevation, distance,
            source_alt2, first_obstruction_angle=0.0;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        /* XXX What to return? */
        return 0.0;
    }
    ReadPath(source,destination,path);

    distance=5280.0*Distance(source,destination);
    source_alt=er+source.alt+GetElevation(source);
    destination_alt=er+destination.alt+GetElevation(destination);
    source_alt2=source_alt*source_alt;

    /* Calculate the cosine of the elevation angle of the
       destination (receiver) as seen by the source (transmitter). */

    cos_xmtr_angle=((source_alt2)+(distance*distance)-(destination_alt*destination_alt))/(2.0*source_alt*distance);

    /* Test all points in between source and destination locations to
       see if the angle to a topographic feature generates a higher
       elevation angle than that produced by the destination.  Begin
       at the source since we're interested in identifying the FIRST
       obstruction along the path between source and destination. */

    for (x=2, block=0; x<path->length && block==0; x++)
    {
        distance=5280.0*path->distance[x];

        test_alt=earthradius+(path->elevation[x]==0.0?path->elevation[x]:path->elevation[x]+clutter);

        cos_test_angle=((source_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*source_alt*distance);

        /* Compare these two angles to determine if
           an obstruction exists.  Since we're comparing
           the cosines of these angles rather than
           the angles themselves, the sense of the
           following "if" statement is reversed from
           what it would be if the angles themselves
           were compared. */

        if (cos_xmtr_angle>=cos_test_angle)
        {
            block=1;
            first_obstruction_angle=((acos(cos_test_angle))/DEG2RAD)-90.0;
        }
    }

    if (block)
        elevation=first_obstruction_angle;

    else
        elevation=((acos(cos_xmtr_angle))/DEG2RAD)-90.0;

    DestroyPath(path);

    return elevation;
}

/* Returns the average terrain calculated in the direction of "azimuth"
 * (degrees) between "start_distance" and "end_distance" (miles) from
 * the source location. If the terrain is all water (non-critical error),
 * -5000.0 is returned. If not enough SDF data has been loaded into
 * memory to complete the survey (critical error), then -9999.0 is
 * returned.
 */
double AverageTerrain(Site source, double azimuthx, double start_distance, double end_distance)
{
    int	c, samples, endpoint;
    double	beta, lat1, lon1, lat2, lon2, num, den, azimuth, terrain=0.0;
    Site destination;

    lat1=source.lat*DEG2RAD;
    lon1=source.lon*DEG2RAD;

    /* Generate a path of elevations between the source
       location and the remote location provided. */

    beta=end_distance/3959.0;

    azimuth=DEG2RAD*azimuthx;

    lat2=asin(sin(lat1)*cos(beta)+cos(azimuth)*sin(beta)*cos(lat1));
    num=cos(beta)-(sin(lat1)*sin(lat2));
    den=cos(lat1)*cos(lat2);

    if (azimuth==0.0 && (beta>HALFPI-lat1))
        lon2=lon1+PI;

    else if (azimuth==HALFPI && (beta>HALFPI+lat1))
        lon2=lon1+PI;

    else if (fabs(num/den)>1.0)
        lon2=lon1;

    else
    {
        if ((PI-azimuth)>=0.0)
            lon2=lon1-arccos(num,den);
        else
            lon2=lon1+arccos(num,den);
    }

    if (lon2<0.0) {
        while (lon2<0.0)
            lon2+=TWOPI;
    }

    if (lon2>TWOPI) {
        while (lon2>TWOPI)
            lon2-=TWOPI;
    }

    lat2=lat2/DEG2RAD;
    lon2=lon2/DEG2RAD;

    destination.lat=lat2;
    destination.lon=lon2;

    /* If SDF data is missing for the endpoint of
       the radial, then the average terrain cannot
       be accurately calculated.  Return -9999.0 */

    if (GetElevation(destination)<-4999.0)
        return (-9999.0);
    else
    {
        Path *path = AllocatePath();
        if (!path) {
            fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
            return (-9999.0);
        }
        ReadPath(source,destination,path);

        endpoint=path->length;

        /* Shrink the length of the radial if the
           outermost portion is not over U.S. land. */

        for (c=endpoint-1; c>=0 && path->elevation[c]==0.0; c--);

        endpoint=c+1;

        for (c=0, samples=0; c<endpoint; c++)
        {
            if (path->distance[c]>=start_distance)
            {
                terrain+=(path->elevation[c]==0.0?path->elevation[c]:path->elevation[c]+clutter);
                samples++;
            }
        }

        if (samples==0)
            terrain=-5000.0;  /* No land */
        else
            terrain=(terrain/(double)samples);

        DestroyPath(path);
        return terrain;
    }
}

/* Returns the antenna's Height Above Average Terrain (HAAT) based on
 * FCC Part 73.313(d). If a critical error occurs, such as a lack of SDF
 * data to complete the survey, -5000.0 is returned.
 */
double haat(Site antenna)
{
    int	azi, c;
    char	error=0;
    double	terrain, avg_terrain, haat, sum=0.0;

    /* Calculate the average terrain between 2 and 10 miles
       from the antenna site at azimuths of 0, 45, 90, 135,
       180, 225, 270, and 315 degrees. */

    for (c=0, azi=0; azi<=315 && error==0; azi+=45)
    {
        terrain=AverageTerrain(antenna, (double)azi, 2.0, 10.0);

        if (terrain<-9998.0)  /* SDF data is missing */
            error=1;

        if (terrain>-4999.0)  /* It's land, not water */
        {
            sum+=terrain;  /* Sum of averages */
            c++;
        }
    }

    if (error)
        return -5000.0;
    else
    {
        avg_terrain=(sum/(double)c);
        if (antenna.amsl_flag)
            haat=antenna.alt-avg_terrain;
        else
            haat=(antenna.alt+GetElevation(antenna))-avg_terrain;

        return haat;
    }
}

/* Places text and marker data in the mask array for illustration on
 * topographic maps generated by SPLAT!. By default, SPLAT! centers text
 * information BELOW the marker, but may move it above, to the left, or
 * to the right of the marker depending on how much room is available on
 * the map, or depending on whether the area is already occupied by
 * another marker or label. If no room or clear space is available on the
 * map to place the marker and its associated text, then the marker and
 * text are not written to the map.
 */
void PlaceMarker(Site location)
{
    int	a, b, c, byte;
    char	ok2print, occupied;
    double	x, y, lat, lon, textx=0.0, texty=0.0, xmin, xmax,
            ymin, ymax, p1, p3, p6, p8, p12, p16, p24, label_length;

    xmin=(double)min_north;
    xmax=(double)max_north;
    ymin=(double)min_west;
    ymax=(double)max_west;
    lat=location.lat;
    lon=location.lon;

    if (lat<xmax && lat>=xmin && (LonDiff(lon,ymax)<=0.0) && (LonDiff(lon,ymin)>=dpp))
    {
        p1=1.0/ppd;
        p3=3.0/ppd;
        p6=6.0/ppd;
        p8=8.0/ppd;
        p12=12.0/ppd;
        p16=16.0/ppd;
        p24=24.0/ppd;

        ok2print=0;
        occupied=0;

        /* Is Marker Position Clear Of Text Or Other Markers? */

        for (a=0, x=lat-p3; (x<=xmax && x>=xmin && a<7); x+=p1, a++)
            for (b=0, y=lon-p3; (LonDiff(y,ymax)<=0.0) && (LonDiff(y,ymin)>=dpp) && b<7; y+=p1, b++)
                occupied|=(GetMask(x,y)&2);

        if (occupied==0)
        {
            /* Determine Where Text Can Be Positioned */

            /* label_length=length in pixels.
               Each character is 8 pixels wide. */

            label_length=p1*(double)(strlen(location.name)<<3);

            if ((LonDiff(lon+label_length,ymax)<=0.0) && (LonDiff(lon-label_length,ymin)>=dpp))
            {
                /* Default: Centered Text */

                texty=lon+label_length/2.0;

                if ((lat-p8)>=p16)
                {
                    /* Position Text Below The Marker */

                    textx=lat-p8;

                    x=textx;
                    y=texty;

                    /* Is This Position Clear Of
                       Text Or Other Markers? */

                    for (a=0, occupied=0; a<16; a++)
                    {
                        for (b=0; b<(int)strlen(location.name); b++)
                            for (c=0; c<8; c++, y-=p1)
                                occupied|=(GetMask(x,y)&2);
                        x-=p1;
                        y=texty;
                    }

                    x=textx;
                    y=texty;

                    if (occupied==0)
                        ok2print=1;
                }

                else
                {
                    /* Position Text Above The Marker */

                    textx=lat+p24;

                    x=textx;
                    y=texty;

                    /* Is This Position Clear Of
                       Text Or Other Markers? */

                    for (a=0, occupied=0; a<16; a++)
                    {
                        for (b=0; b<(int)strlen(location.name); b++)
                            for (c=0; c<8; c++, y-=p1)
                                occupied|=(GetMask(x,y)&2);
                        x-=p1;
                        y=texty;
                    }

                    x=textx;
                    y=texty;

                    if (occupied==0)
                        ok2print=1;
                }
            }

            if (ok2print==0)
            {
                if (LonDiff(lon-label_length,ymin)>=dpp)
                {
                    /* Position Text To The
                       Right Of The Marker */

                    textx=lat+p6;
                    texty=lon-p12;

                    x=textx;
                    y=texty;

                    /* Is This Position Clear Of
                       Text Or Other Markers? */

                    for (a=0, occupied=0; a<16; a++)
                    {
                        for (b=0; b<(int)strlen(location.name); b++)
                            for (c=0; c<8; c++, y-=p1)
                                occupied|=(GetMask(x,y)&2);
                        x-=p1;
                        y=texty;
                    }

                    x=textx;
                    y=texty;

                    if (occupied==0)
                        ok2print=1;
                }

                else
                {
                    /* Position Text To The
                       Left Of The Marker */

                    textx=lat+p6;
                    texty=lon+p8+(label_length);

                    x=textx;
                    y=texty;

                    /* Is This Position Clear Of
                       Text Or Other Markers? */

                    for (a=0, occupied=0; a<16; a++)
                    {
                        for (b=0; b<(int)strlen(location.name); b++)
                            for (c=0; c<8; c++, y-=p1)
                                occupied|=(GetMask(x,y)&2);
                        x-=p1;
                        y=texty;
                    }

                    x=textx;
                    y=texty;

                    if (occupied==0)
                        ok2print=1;
                }
            }

            /* textx and texty contain the latitude and longitude
               coordinates that describe the placement of the text
               on the map. */

            if (ok2print)
            {
                /* Draw Text */

                x=textx;
                y=texty;

                for (a=0; a<16; a++)
                {
                    for (b=0; b<(int)strlen(location.name); b++)
                    {
                        byte=fontdata[16*(location.name[b])+a];

                        for (c=128; c>0; c=c>>1, y-=p1)
                            if (byte&c)
                                OrMask(x,y,2);
                    }

                    x-=p1;
                    y=texty;
                }

                /* Draw Square Marker Centered
                   On Location Specified */

                for (a=0, x=lat-p3; (x<=xmax && x>=xmin && a<7); x+=p1, a++)
                    for (b=0, y=lon-p3; (LonDiff(y,ymax)<=0.0) && (LonDiff(y,ymin)>=dpp) && b<7; y+=p1, b++)
                        OrMask(x,y,2);
            }
        }
    }
}

/* This function takes numeric input in the form of a character
 * string, and returns an equivalent bearing in degrees as a
 * decimal number (double).  The input may either be expressed
 * in decimal format (40.139722) or degree, minute, second
 * format (40 08 23).  This function also safely handles
 * extra spaces found either leading, trailing, or
 * embedded within the numbers expressed in the
 * input string.  Decimal seconds are permitted.
 */
double ReadBearing(char *input)
{
    double	seconds, bearing=0.0;
    char	buf[20];
    int	a, b, length, degrees, minutes;

    /* Copy "input" to "buf", and ignore any extra
       spaces that might be present in the process. */

    buf[0]=0;
    length= (int)strlen(input);

    for (a=0, b=0; a<length && a<18; a++)
    {
        if ((input[a]!=32 && input[a]!='\n') || (input[a]==32 && input[a+1]!=32 && input[a+1]!='\n' && b!=0))
        {
            buf[b]=input[a];
            b++;
        }	 
    }

    buf[b]=0;

    /* Count number of spaces in the clean string. */

    length= (int)strlen(buf);

    for (a=0, b=0; a<length; a++)
        if (buf[a]==32)
            b++;

    if (b==0)  /* Decimal Format (40.139722) */
        sscanf(buf,"%lf",&bearing);

    if (b==2)  /* Degree, Minute, Second Format (40 08 23.xx) */
    {
        sscanf(buf,"%d %d %lf",&degrees, &minutes, &seconds);

        bearing=fabs((double)degrees);
        bearing+=fabs(((double)minutes)/60.0);
        bearing+=fabs(seconds/3600.0);

        if ((degrees<0) || (minutes<0) || (seconds<0.0))
            bearing=-bearing;
    }

    /* Anything else returns a 0.0 */

    if (bearing>360.0 || bearing<-360.0)
        bearing=0.0;

    return bearing;
}

/* Reads SPLAT! .qth (site location) files. The latitude and
 * longitude may be expressed either in decimal degrees, or in
 * degree, minute, second format.  
 * 
 * Antenna height is assumed to be expressed in feet above
 * ground level (AGL), unless followed by the letter 'M',
 * or 'm', or by the word "meters", "Meters", or "METERS",
 * in which case meters is assumed, and is converted to feet
 * before exiting.  Additionally, if the height given is
 * followed by "amsl" or "AMSL", it is recognized as being
 * the height above mean sea level rather than the height
 * above ground (AGL), and the amsl_flag is set accordingly.
 */
Site LoadQTH(char *filename)
{
    int	x;
    char	buf[MAX_LINE_LEN], qthfile[MAX_PATH_LEN];
    Site    tempsite;
    FILE	*fd=NULL;

    /* XXX FIX ME */
    x= (int)strlen(filename);
    strncpy(qthfile, filename, MAX_PATH_LEN-1);

    if (qthfile[x-3]!='q' || qthfile[x-2]!='t' || qthfile[x-1]!='h')
    {
        if (x>MAX_PATH_LEN-4)
            qthfile[MAX_PATH_LEN-4]=0;

        strcat(qthfile,".qth");
    }


    tempsite.lat=91.0;
    tempsite.lon=361.0;
    tempsite.alt=0.0;
    tempsite.name[0]=0;
    tempsite.filename[0]=0;

    fd=fopen(qthfile,"r");

    if (fd!=NULL)
    {
        /* Site Name */
        fgets(buf,MAX_LINE_LEN,fd);

        /* Strip <CR> and/or <LF> from end of site name */

        for (x=0; buf[x]!=13 && buf[x]!=10 && buf[x]!=0 && x<(MAX_SITE_NAME_LEN-1); tempsite.name[x]=buf[x], x++);

        tempsite.name[x]=0;

        /* Site Latitude */
        fgets(buf,MAX_LINE_LEN,fd);
        tempsite.lat=ReadBearing(buf);

        /* Site Longitude */
        fgets(buf,MAX_LINE_LEN,fd);
        tempsite.lon=ReadBearing(buf);

        if (tempsite.lon<0.0)
            tempsite.lon+=360.0;

        /* Antenna Height */
        fgets(buf,MAX_LINE_LEN,fd);
        fclose(fd);

        /* Remove <CR> and/or <LF> from antenna height
           buf and convert any related text to uppercase */

        for (x=0; buf[x]!=13 && buf[x]!=10 && buf[x]!=0; buf[x]=toupper(buf[x]), x++);

        buf[x]=0;

        /* Read antenna height from buf */

        sscanf(buf,"%lf",&tempsite.alt);

        if ((buf[x-1]=='M') || strstr(buf,"M ") || strstr(buf,"METERS"))
            tempsite.alt*=3.28084;

        if (strstr(buf,"AMSL"))
            tempsite.amsl_flag=1;
        else
            tempsite.amsl_flag=0;

        for (x=0; x<(MAX_PATH_LEN-1) && qthfile[x]!=0; x++)
            tempsite.filename[x]=qthfile[x];

        tempsite.filename[x]=0;
    }

    return tempsite;
}

/* Reads and processes antenna pattern (.az and .el) files that correspond
 * in name to previously loaded SPLAT! .lrp files.
 */
void LoadPAT(char *filename)
{
    int	a, b, w, x, y, z, last_index, next_index, span;
    char buf[MAX_LINE_LEN], *pointer=NULL;
    char *azfile, *elfile;
    float	az, xx, amplitude, rotation, valid1, valid2,
            delta, azimuth[361], azimuth_pattern[361],
            mechanical_tilt = 0.0, tilt_azimuth, sum;
    double elevation;
    /* elevation_pattern[361][1001]; */ /* 2.7MB. This blows the stack right out of the water. Allocate on heap. */
    double tilt, tilt_increment, slant_angle[361];
    FILE	*fd=NULL;

    float *el_pattern = (float*)malloc(sizeof(float) * 10001);
    unsigned char *read_count = (unsigned char*)malloc(sizeof(unsigned char*) * 10001);
    double **elevation_pattern = (double**)malloc(sizeof(double*) * 361);
    for (x = 0; x < 361; ++x) {
        elevation_pattern[x] = (double*)malloc(sizeof(double) * 1001);
    }

    azfile = copyFilename(filename, "az");
    elfile = copyFilename(filename, "el");

    rotation=0.0;

    got_azimuth_pattern=false;
    got_elevation_pattern=false;

    /* Load .az antenna pattern file */
    fd=fopen(azfile,"r");

    if (fd!=NULL)
    {
        printf("Loading %s\n", azfile);

        /* Clear azimuth pattern array */

        for (x=0; x<=360; x++)
        {
            azimuth[x]=0.0;
            read_count[x]=0;
        }


        /* Read azimuth pattern rotation
           in degrees measured clockwise
           from true North. */

        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(buf,"%f",&rotation);


        /* Read azimuth (degrees) and corresponding
           normalized field radiation pattern amplitude
           (0.0 to 1.0) until EOF is reached. 

           NOTES
           1) If a radial has more than one entry, the azimuth[x]
           for each will be added. This is averaged out later.
           2) It will interpolate to fill in any missing radials, but will
           not do this past the 360'th radial.
           3) There is no check to make sure there are at least two
           radials entered (which is necessary for the averaging).
         */

        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(buf,"%f %f",&az, &amplitude);

        do
        {
            x=(int)rintf(az);

            if (x>=0 && x<=360 && fd!=NULL)
            {
                azimuth[x]+=amplitude;
                read_count[x]++;
            }

            fgets(buf,MAX_LINE_LEN,fd);
            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            sscanf(buf,"%f %f",&az, &amplitude);

        } while (feof(fd)==0);

        fclose(fd);


        /* Handle 0=360 degree ambiguity */

        if ((read_count[0]==0) && (read_count[360]!=0))
        {
            read_count[0]=read_count[360];
            azimuth[0]=azimuth[360];
        }

        if ((read_count[0]!=0) && (read_count[360]==0))
        {
            read_count[360]=read_count[0];
            azimuth[360]=azimuth[0];
        }

        /* Average pattern values in case more than
           one was read for each degree of azimuth. */
        /* Also, check count of read values */

        for (sum=0, x=0; x<=360; x++)
        {
            if (read_count[x]>1) {
                azimuth[x]/=(float)read_count[x];
            }
            if (read_count[x]!=0) {
                sum++;
            }
        }

        if (sum==0) {
            printf("*** WARNING: There were no radials in the file. This probably isn't correct.\n");
        } else if (sum==1) {
            printf("*** WARNING: There was only one radial in the file. No interpolation can be done.\n");
        }

        /* Interpolate missing azimuths
           to completely fill the array */

        last_index=-1;
        next_index=-1;

        for (x=0; x<=360; x++)
        {
            if (read_count[x]!=0)
            {
                if (last_index==-1)
                    last_index=x;
                else
                    next_index=x;
            }

            if (last_index!=-1 && next_index!=-1)
            {
                valid1=azimuth[last_index];
                valid2=azimuth[next_index];

                span=next_index-last_index;
                delta=(valid2-valid1)/(float)span;

                for (y=last_index+1; y<next_index; y++)
                    azimuth[y]=azimuth[y-1]+delta;

                last_index=y;
                next_index=-1;
            }
        }


        /* Perform azimuth pattern rotation
           and load azimuth_pattern[361] with
           azimuth pattern data in its final form. */

        for (x=0; x<360; x++)
        {
            y=x+(int)rintf(rotation);

            if (y>=360)
                y-=360;

            azimuth_pattern[y]=azimuth[x];
        }

        azimuth_pattern[360]=azimuth_pattern[0];

        got_azimuth_pattern=true;
    }
    free(azfile);

    /* Read and process .el file */

    fd=fopen(elfile,"r");

    if (fd!=NULL)
    {
        printf("Loading %s\n", elfile);

        for (x=0; x<=10000; x++)
        {
            el_pattern[x]=0.0;
            read_count[x]=0;
        }

        /* Read mechanical tilt (degrees) and
           tilt azimuth in degrees measured
           clockwise from true North. */  

        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(buf,"%f %f",&mechanical_tilt, &tilt_azimuth);

        /* Read elevation (degrees) and corresponding
           normalized field radiation pattern amplitude
           (0.0 to 1.0) until EOF is reached.

           NOTES
           1) If a radial has more than one entry, the azimuth[x]
           for each will be added. This is averaged out later.
           2) It will interpolate to fill in any missing radials, but will
           not do this past the 360'th radial.
           3) There is no check to make sure there are at least two
           radials entered (which is necessary for the averaging).
         */

        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(buf,"%lf %f", &elevation, &amplitude);

        while (feof(fd)==0)
        {
            /* Read in normalized radiated field values
               for every 0.01 degrees of elevation between
               -10.0 and +90.0 degrees */

            x=(int)rint(100.0*(elevation+10.0));

            if (x>=0 && x<=10000)
            {
                el_pattern[x]+=amplitude;
                read_count[x]++;
            }

            fgets(buf,254,fd);
            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            sscanf(buf,"%lf %f", &elevation, &amplitude);
        }

        fclose(fd);

        /* Average the field values in case more than
           one was read for each 0.01 degrees of elevation. */
        /* Also, check count of read values */

        for (sum=0,x=0; x<=10000; x++)
        {
            if (read_count[x]>1) {
                el_pattern[x]/=(float)read_count[x];
            }
            if (read_count[x]!=0) {
                sum++;
            }
        }

        if (sum==0) {
            printf("*** WARNING: There were no elevation angles in the file. This probably isn't correct.\n");
        } else if (sum==1) {
            printf("*** WARNING: There was only one elevation angle in the file. No interpolation was done.\n");
        }


        /* Interpolate between missing elevations (if
           any) to completely fill the array and provide
           radiated field values for every 0.01 degrees of
           elevation. */

        last_index=-1;
        next_index=-1;

        for (x=0; x<=10000; x++)
        {
            if (read_count[x]!=0)
            {
                if (last_index==-1)
                    last_index=x;
                else
                    next_index=x;
            }

            if (last_index!=-1 && next_index!=-1)
            {
                valid1=el_pattern[last_index];
                valid2=el_pattern[next_index];

                span=next_index-last_index;
                delta=(valid2-valid1)/(float)span;

                for (y=last_index+1; y<next_index; y++)
                    el_pattern[y]=el_pattern[y-1]+delta;

                last_index=y;
                next_index=-1;
            }
        }

        /* Fill slant_angle[] array with offset angles based
           on the antenna's mechanical beam tilt (if any)
           and tilt direction (azimuth). */

        if (mechanical_tilt==0.0)
        {
            for (x=0; x<=360; x++)
                slant_angle[x]=0.0;
        }

        else
        {
            tilt_increment=mechanical_tilt/90.0;

            for (x=0; x<=360; x++)
            {
                xx=(float)x;
                y=(int)rintf(tilt_azimuth+xx);

                while (y>=360)
                    y-=360;

                while (y<0)
                    y+=360;

                if (x<=180)
                    slant_angle[y]=-(tilt_increment*(90.0-xx));

                if (x>180)
                    slant_angle[y]=-(tilt_increment*(xx-270.0));
            }
        }

        slant_angle[360]=slant_angle[0];   /* 360 degree wrap-around */

        for (w=0; w<=360; w++)
        {
            tilt=slant_angle[w];

            /** Convert tilt angle to
              an array index offset **/

            y=(int)rint(100.0*tilt);

            /* Copy shifted el_pattern[10001] field
               values into elevation_pattern[361][1001]
               at the corresponding azimuth, downsampling
               (averaging) along the way in chunks of 10. */

            for (x=y, z=0; z<=1000; x+=10, z++)
            {
                for (sum=0.0, a=0; a<10; a++)
                {
                    b=a+x;

                    if (b>=0 && b<=10000)
                        sum+=el_pattern[b];
                    if (b<0)
                        sum+=el_pattern[0];
                    if (b>10000)
                        sum+=el_pattern[10000];
                }

                elevation_pattern[w][z]=sum/10.0;
            }
        }

        got_elevation_pattern=true;
    }
    free(elfile);

    for (x=0; x<=360; x++)
    {
        for (y=0; y<=1000; y++)
        {
            if (got_elevation_pattern)
                elevation=elevation_pattern[x][y];
            else
                elevation=1.0;

            if (got_azimuth_pattern)
                az=azimuth_pattern[x];
            else
                az=1.0;

            LR.antenna_pattern[x][y]=az*elevation;
        }
    }

    for (x = 0; x < 361; ++x) {
        free(elevation_pattern[x]);
    }
    free(elevation_pattern);
    free(el_pattern);
    free(read_count);
}

/* Returns at most one less than 'length' number of characters
 * from * a bz2 compressed file whose file descriptor is pointed
 * to by *bzfp. In operation, a buffer is filled with uncompressed
 * data (size = BZBUFFER), which is then parsed and doled out
 * as NULL-terminated character bufs every time this function is
 * invoked. A NULL buf indicates an EOF or error condition.
 *
 * Uses internal static buffers. NOT THREADSAFE.
 */
char *BZfgets(BZFILE *bzfp, unsigned length)
{
    static int x, y, nBuf;
    static char buffer[BZBUFFER+1], output[BZBUFFER+1];
    char done=0;

    if (!opened && bzerror==BZ_OK)
    {
        /* First time through.  Initialize everything! */

        x=0;
        y=0;
        nBuf=0;
        opened=true;
        output[0]=0;
    }

    do
    {
        if (x==nBuf && bzerror!=BZ_STREAM_END && bzerror==BZ_OK && opened)
        {
            /* Uncompress data into a static buffer */

            nBuf=BZ2_bzRead(&bzerror, bzfp, buffer, BZBUFFER);
            buffer[nBuf]=0;
            x=0;
        }

        /* Build a buf from buffer contents */

        output[y]=buffer[x];

        if (output[y]=='\n' || output[y]==0 || y==(int)length-1)
        {
            output[y+1]=0;
            done=1;
            y=0;
        }

        else
            y++;
        x++;

    } while (done==0);

    if (output[0]==0)
        opened=false;

    return (output);
}

/* Reads SPLAT Data Files (.sdf) containing digital elevation model data
 * into memory. It loads them into DEM structs, which then get appended
 * to the global dem array.
 */
int LoadSDF(int minlat, int maxlat, int minlon, int maxlon, bool hidef)
{
    int	x, y, data;
    char line[64];
    char sdf_file[MAX_PATH_LEN*2];
    SDFCompressType compressType = SDF_COMPRESSTYPE_NONE;
    DEM *dem;
    FILE *fp = NULL;
    BZFILE *bzfp = NULL;

    /* this sets both the kinds of formats we understand and the priority */
    SDFCompressFormat formats[] = {
        { SDF_COMPRESSTYPE_NONE, "sdf" },
        { SDF_COMPRESSTYPE_BZIP2, "sdf.bz2" },
    };
    const int known_formats = sizeof(formats)/sizeof(SDFCompressFormat);

    /* Is it already in memory? */

    dem = FindDEM_Explicit(minlat, minlon, maxlat, maxlon);
    if (dem != NULL) {
        return 0;
    }

    dem = AllocateDEM();
    if (dem == NULL) {
        return -1;
    }

    /* Rearranged priority! Old way was:
     *	1) look in current workdir first
     *	2) look in directory specified by -d argument or .splat_path files
     *
     * New way:
     *	1) Look in directory specified by -d argument
     *  2) Look in local directory
     *  3) Look in directory named by $HOME/.splat_path
     *
     * My reasoning is that if someone actively specifies something on the
     * command line, they probably mean it to override any other settings.
     *
     */

    /* check -d directory */
    if (sdf_path[0]!=0) {
        for (int i=0; i<known_formats; ++i) {
            snprintf(sdf_file, MAX_PATH_LEN*2, "%s%d_%d_%d_%d%s.%s", sdf_path,
                     minlat, maxlat, minlon, maxlon, (hidef?"-hd":""), formats[i].suffix);
            if ((fp=fopen(sdf_file,"rb"))!= NULL) {
                compressType=formats[i].type;
                break;
            }
#ifndef _WIN32
            snprintf(sdf_file, MAX_PATH_LEN*2, "%s%d:%d:%d:%d%s.%s", sdf_path,
                     minlat, maxlat, minlon, maxlon, (hidef?"-hd":""), formats[i].suffix);
            if ((fp=fopen(sdf_file,"rb"))!= NULL) {
                compressType=formats[i].type;
                break;
            }
#endif
        }
    }

    /* check local directory */
    if (!fp) {
        for (int i=0; i<known_formats; ++i) {
            snprintf(sdf_file, MAX_PATH_LEN*2, "%d_%d_%d_%d%s.%s",
                     minlat, maxlat, minlon, maxlon, (hidef?"-hd":""), formats[i].suffix);
            if ((fp=fopen(sdf_file,"rb"))!= NULL) {
                compressType=formats[i].type;
                break;
            }
#ifndef _WIN32
            snprintf(sdf_file, MAX_PATH_LEN*2, "%d:%d:%d:%d%s.%s",
                     minlat, maxlat, minlon, maxlon, (hidef?"-hd":""), formats[i].suffix);
            if ((fp=fopen(sdf_file,"rb"))!= NULL) {
                compressType=formats[i].type;
                break;
            }
#endif
        }
    }

    /* check $HOME/.splat_path directory */
    if (!fp && home_sdf_path[0]!=0) {
        for (int i=0; i<known_formats; ++i) {
            snprintf(sdf_file, MAX_PATH_LEN*2, "%s%d_%d_%d_%d%s.%s", home_sdf_path,
                     minlat, maxlat, minlon, maxlon, (hidef?"-hd":""), formats[i].suffix);
            if ((fp=fopen(sdf_file,"rb"))!= NULL) {
                compressType=formats[i].type;
                break;
            }
#ifndef _WIN32
            snprintf(sdf_file, MAX_PATH_LEN*2, "%s%d:%d:%d:%d%s.%s", home_sdf_path,
                     minlat, maxlat, minlon, maxlon, (hidef?"-hd":""), formats[i].suffix);
            if ((fp=fopen(sdf_file,"rb"))!= NULL) {
                compressType=formats[i].type;
                break;
            }
#endif
        }
    }

    if (fp) {
        fprintf(stdout,"Loading \"%s\"...", sdf_file);
        fflush(stdout);

        if (compressType == SDF_COMPRESSTYPE_BZIP2) {
            bzfp=BZ2_bzReadOpen(&bzerror,fp,0,0,NULL,0);
        }

        switch (compressType) {
            case SDF_COMPRESSTYPE_BZIP2:
                dem->max_west = atoi(BZfgets(bzfp, 64));
                dem->min_north = atoi(BZfgets(bzfp, 64));
                dem->min_west = atoi(BZfgets(bzfp, 64));
                dem->max_north = atoi(BZfgets(bzfp, 64));
                break;
            case SDF_COMPRESSTYPE_NONE:
            default:
                dem->max_west = atoi(fgets(line, 64, fp));
                dem->min_north = atoi(fgets(line, 64, fp));
                dem->min_west = atoi(fgets(line, 64, fp));
                dem->max_north = atoi(fgets(line, 64, fp));
        }

        for (x=0; x<ippd; x++) {
            for (y=0; y<ippd; y++) {
                switch (compressType) {
                    case SDF_COMPRESSTYPE_BZIP2:
                        data = atoi(BZfgets(bzfp, 64));
                        break;
                    case SDF_COMPRESSTYPE_NONE:
                    default:
                        data = atoi(fgets(line, 64, fp));
                }

                dem->data[x][y]=data;
                dem->signal[x][y]=0;
                dem->mask[x][y]=0;

                if (data>dem->max_el)
                    dem->max_el=data;

                if (data<dem->min_el)
                    dem->min_el=data;
            }
        }

        if (compressType == SDF_COMPRESSTYPE_BZIP2) {
            BZ2_bzReadClose(&bzerror,bzfp);
        }
        fclose(fp);

    } else {
        fprintf(stdout,"Can't find file for region %d_%d_%d_%d. Assuming as sea-level.\n",
                minlat, maxlat, minlon, maxlon);
        fflush(stdout);

        dem->max_west=maxlon;
        dem->min_north=minlat;
        dem->min_west=minlon;
        dem->max_north=maxlat;
        dem->max_el=0;
        dem->min_el=0;

        for (x=0; x<ippd; x++) {
            for (y=0; y<ippd; y++) {
                dem->data[x][y]=0;
                dem->signal[x][y]=0;
                dem->mask[x][y]=0;
            }
        }

    }

    if (dem->min_el<min_elevation)
        min_elevation=dem->min_el;

    if (dem->max_el>max_elevation)
        max_elevation=dem->max_el;

    if (max_north==-90)
        max_north=dem->max_north;

    else if (dem->max_north>max_north)
        max_north=dem->max_north;

    if (min_north==90)
        min_north=dem->min_north;

    else if (dem->min_north<min_north)
        min_north=dem->min_north;

    if (max_west==-1) {
        max_west=dem->max_west;
    } else {
        if (abs(dem->max_west-max_west)<180) {
            if (dem->max_west>max_west)
                max_west=dem->max_west;
        } else {
            if (dem->max_west<max_west)
                max_west=dem->max_west;
        }
    }

    if (min_west==360) {
        min_west=dem->min_west;
    } else {
        if (abs(dem->min_west-min_west)<180) {
            if (dem->min_west<min_west)
                min_west=dem->min_west;
        } else {
            if (dem->min_west>min_west)
                min_west=dem->min_west;
        }
    }

    AppendDEM(dem); /* append to global array */

    fprintf(stdout," Done!\n");
    fflush(stdout);

    return 1;
}

/* This function reads SPLAT! city/site files, and plots the locations
 * and names of the cities and site locations read on topographic maps
 * generated by SPLAT!
 */
void LoadCities(char *filename)
{
    int	x, y, z;
    char	input[MAX_LINE_LEN], str[3][MAX_LINE_LEN];
    Site	city_site;
    FILE	*fd=NULL;

    fd=fopen(filename,"r");

    if (fd!=NULL)
    {
        fgets(input,MAX_LINE_LEN,fd);

        fprintf(stdout,"\nReading \"%s\"... ",filename);
        fflush(stdout);

        while (fd!=NULL && feof(fd)==0)
        {
            /* Parse line for name, latitude, and longitude */

            for (x=0, y=0, z=0; x<MAX_LINE_LEN && input[x]!=0 && z<3; x++)
            {
                if (input[x]!=',' && y<MAX_LINE_LEN)
                {
                    str[z][y]=input[x];
                    y++;
                }
                else
                {
                    str[z][y]='\0';
                    z++;
                    y=0;
                }
            }

            strncpy(city_site.name,str[0],MAX_SITE_NAME_LEN-1);
            city_site.name[MAX_SITE_NAME_LEN-1] = '\0';
            city_site.lat=ReadBearing(str[1]);
            city_site.lon=ReadBearing(str[2]);
            city_site.alt=0.0;

            if (city_site.lon<0.0)
                city_site.lon+=360.0;

            PlaceMarker(city_site);

            fgets(input,MAX_LINE_LEN,fd);
        }

        fclose(fd);
        fprintf(stdout,"Done!");
        fflush(stdout);
    }

    else
        fprintf(stderr,"\n*** ERROR: \"%s\": not found!",filename);
}

#ifndef _WIN32
/* Reads a file containing User-Defined Terrain features for their
 * addition to the digital elevation model data used by SPLAT!. 
 * Elevations in the UDT file are evaluated and then copied into a
 * temporary file under /tmp. Then the contents of the temp file are
 * scanned, and if found to be unique, are added to the ground
 * elevations described by the digital elevation data already loaded
 * into memory.
 *
 * Doesn't currently work on Windows, primiarly because of the use
 * of mkstemp.
 */
void LoadUDT(char *filename)
{
    int	i, x, y, z, ypix, xpix, tempxpix, tempypix, fd=0;
    char	input[MAX_LINE_LEN], str[3][MAX_LINE_LEN], tempname[15], *pointer=NULL;
    double	latitude, longitude, height, tempheight;
    FILE	*fd1=NULL, *fd2=NULL;

    strcpy(tempname,"/tmp/XXXXXX\0");

    fd1=fopen(filename,"r");

    if (fd1!=NULL)
    {
        fd=mkstemp(tempname);
        fd2=fopen(tempname,"w");

        fgets(input,MAX_LINE_LEN,fd1);

        pointer=strchr(input,';');

        if (pointer!=NULL)
            *pointer=0;

        fprintf(stdout,"\nReading \"%s\"... ",filename);
        fflush(stdout);

        while (feof(fd1)==0)
        {
            /* Parse line for latitude, longitude, height */

            for (x=0, y=0, z=0; x<MAX_LINE_LEN && input[x]!=0 && z<3; x++)
            {
                if (input[x]!=',' && y<MAX_LINE_LEN)
                {
                    str[z][y]=input[x];
                    y++;
                }

                else
                {
                    str[z][y]='\0';
                    z++;
                    y=0;
                }
            }

            latitude=ReadBearing(str[0]);
            longitude=ReadBearing(str[1]);

            if (longitude<0.0)
                longitude+=360.0;

            /* Remove <CR> and/or <LF> from antenna height buf */

            for (i=0; str[2][i]!=13 && str[2][i]!=10 && str[2][i]!=0; i++);

            str[2][i]=0;

            /* The terrain feature may be expressed in either
               feet or meters.  If the letter 'M' or 'm' is
               discovered in the buf, then this is an
               indication that the value given is expressed
               in meters.  Otherwise the height is interpreted
               as being expressed in feet.  */

            for (i=0; str[2][i]!='M' && str[2][i]!='m' && str[2][i]!=0 && i<MAX_LINE_LEN; i++);

            if (str[2][i]=='M' || str[2][i]=='m')
            {
                str[2][i]=0;
                height=rint(atof(str[2]));
            }

            else
            {
                str[2][i]=0;
                height=rint(METERS_PER_FOOT*atof(str[2]));
            }

            if (height>0.0)
                fprintf(fd2,"%d, %d, %f\n",(int)rint(latitude/dpp), (int)rint(longitude/dpp), height);

            fgets(input,MAX_LINE_LEN,fd1);

            pointer=strchr(input,';');

            if (pointer!=NULL)
                *pointer=0;
        }

        fclose(fd1);
        fclose(fd2);
        close(fd);

        fprintf(stdout,"Done!");
        fflush(stdout);

        fd1=fopen(tempname,"r");
        fd2=fopen(tempname,"r");

        y=0;

        fscanf(fd1,"%d, %d, %lf", &xpix, &ypix, &height);

        do
        {
            x=0;
            z=0;

            fscanf(fd2,"%d, %d, %lf", &tempxpix, &tempypix, &tempheight);

            do
            {
                if (x>y && xpix==tempxpix && ypix==tempypix)
                {
                    z=1;  /* Dupe! */

                    if (tempheight>height)
                        height=tempheight;
                }

                else
                {
                    fscanf(fd2,"%d, %d, %lf", &tempxpix, &tempypix, &tempheight);
                    x++;
                }

            } while (feof(fd2)==0 && z==0);

            if (z==0)  /* No duplicate found */
                AddElevation(xpix*dpp, ypix*dpp, height);

            fscanf(fd1,"%d, %d, %lf", &xpix, &ypix, &height);
            y++;

            rewind(fd2);

        } while (feof(fd1)==0);

        fclose(fd1);
        fclose(fd2);
        unlink(tempname);
    }

    else
        fprintf(stderr,"\n*** ERROR: \"%s\": not found!",filename);

    fprintf(stdout,"\n");
}
#endif

/* Reads Cartographic Boundary Files available from the U.S. Census Bureau,
 * and plots the data contained in those files on the map images generated
 * by SPLAT!. Such files contain the coordinates that describe the boundaries
 * of cities, counties, and states.
 */
void LoadBoundaries(char *filename)
{
    int	x;
    double	lat0, lon0, lat1, lon1;
    char	buf[MAX_LINE_LEN];
    Site    source, destination;
    FILE	*fd=NULL;

    fd=fopen(filename,"r");

    if (fd!=NULL)
    {
        fgets(buf,MAX_LINE_LEN,fd);

        fprintf(stdout,"\nReading \"%s\"... ",filename);
        fflush(stdout);

        do
        {
            fgets(buf,MAX_LINE_LEN,fd);
            sscanf(buf,"%lf %lf", &lon0, &lat0);

            fgets(buf,MAX_LINE_LEN,fd);
            do
            {
                sscanf(buf,"%lf %lf", &lon1, &lat1);

                source.lat=lat0;
                source.lon=(lon0>0.0 ? 360.0-lon0 : -lon0);
                destination.lat=lat1;
                destination.lon=(lon1>0.0 ? 360.0-lon1 : -lon1);

                Path *path = AllocatePath();
                if (!path) {
                    fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
                    return;
                }
                ReadPath(source,destination,path);

                for (x=0; x<path->length; x++)
                    OrMask(path->lat[x],path->lon[x],4);

                DestroyPath(path);

                lat0=lat1;
                lon0=lon1;

                fgets(buf,78,fd);

            } while (strncmp(buf,"END",3)!=0 && feof(fd)==0);

            fgets(buf,MAX_LINE_LEN,fd);

        } while (strncmp(buf,"END",3)!=0 && feof(fd)==0);

        fclose(fd);

        fprintf(stdout,"Done!");
        fflush(stdout);
    }

    else
        fprintf(stderr,"\n*** ERROR: \"%s\": not found!",filename);
}

/* Reads ITM parameter data for the transmitter site. The file name
 * is the same as the txsite, except the filename extension is .lrp.
 * If the needed file is not found, then the file "splat.lrp" is read
 * from the current working directory. Failure to load this file under
 * a forced_read condition will result in the default parameters
 * hard-coded into this function to be used and written to "splat.lrp".
 */
char LoadLRP(Site txsite, char forced_read)
{
    double	din;
    char buf[MAX_LINE_LEN], *pointer=NULL, return_value=0;
    char *filename;
    int	iin, ok=0;
    FILE	*fd=NULL, *outfile=NULL;

    /* Default parameters */

    LR.eps_dielect=0.0;
    LR.sgm_conductivity=0.0;
    LR.eno_ns_surfref=0.0;
    LR.frq_mhz=0.0;
    LR.radio_climate=0;
    LR.pol=0;
    LR.conf=0.0;
    LR.rel=0.0;
    LR.erp=0.0;

    for (int i = 0; i < ANTENNA_RADIALS; ++i) {
        for (int j = 0; j < ANTENNA_ANGLES; ++j) {
            LR.antenna_pattern[i][j] = NO_ANTENNA_DATA;
        }
    }

    /* Generate .lrp filename from txsite filename. */

    filename = copyFilename(txsite.filename, "lrp");

    printf("Loading %s\n", filename);

    fd=fopen(filename,"r");

    if (fd==NULL)
    {
        free(filename);

        /* Load default "splat.lrp" file */
        filename = strdup("splat.lrp");
        fd=fopen(filename,"r");
    }

    if (fd!=NULL)
    {
        fgets(buf,MAX_LINE_LEN,fd);

        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        ok=sscanf(buf,"%lf", &din);

        if (ok)
        {
            LR.eps_dielect=din;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%lf", &din);
        }

        if (ok)
        {
            LR.sgm_conductivity=din;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%lf", &din);
        }

        if (ok)
        {
            LR.eno_ns_surfref=din;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%lf", &din);
        }

        if (ok)
        {
            LR.frq_mhz=din;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%d", &iin);
        }

        if (ok)
        {
            LR.radio_climate=iin;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%d", &iin);
        }

        if (ok)
        {
            LR.pol=iin;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%lf", &din);
        }

        if (ok)
        {
            LR.conf=din;

            fgets(buf,MAX_LINE_LEN,fd);

            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%lf", &din);
        }

        if (ok)
        {
            LR.rel=din;
            din=0.0;
            return_value=1;

            if (fgets(buf,MAX_LINE_LEN,fd)!=NULL)
            {
                pointer=strchr(buf,';');

                if (pointer!=NULL)
                    *pointer=0;

                if (sscanf(buf,"%lf", &din))
                    LR.erp=din;

                /* ERP in SPLAT! is referenced to 1 Watt
                   into a dipole (0 dBd).  If ERP is
                   expressed in dBm (referenced to a
                   0 dBi radiator), convert dBm in EIRP
                   to ERP.  */

                if ((strstr(buf, "dBm")!=NULL) || (strstr(buf,"dbm")!=NULL))
                    LR.erp=(pow(10.0,(LR.erp-32.14)/10.0));
            }
        }

        fclose(fd);

        if (forced_erp!=-1.0)
            LR.erp=forced_erp;

        if (forced_freq>=20.0 && forced_freq<=20000.0)
            LR.frq_mhz=forced_freq;

        if (ok)
            LoadPAT(filename);
    }

    if (fd==NULL && forced_read)
    {
        /* Assign some default parameters
           for use in this run. */

        LR.eps_dielect=15.0;
        LR.sgm_conductivity=0.005;
        LR.eno_ns_surfref=301.0;
        LR.frq_mhz=300.0;
        LR.radio_climate=5;
        LR.pol=0;
        LR.conf=0.50;
        LR.rel=0.50;
        LR.erp=0.0;

        /* Write them to a "splat.lrp" file. */

        outfile=fopen("splat.lrp","w");

        fprintf(outfile,"%.3f\t; Earth Dielectric Constant (Relative permittivity)\n",LR.eps_dielect);
        fprintf(outfile,"%.3f\t; Earth Conductivity (Siemens per meter)\n", LR.sgm_conductivity);
        fprintf(outfile,"%.3f\t; Atmospheric Bending Constant (N-Units)\n",LR.eno_ns_surfref);
        fprintf(outfile,"%.3f\t; Frequency in MHz (20 MHz to 20 GHz)\n", LR.frq_mhz);
        fprintf(outfile,"%d\t; Radio Climate\n",LR.radio_climate);
        fprintf(outfile,"%d\t; Polarization (0 = Horizontal, 1 = Vertical)\n", LR.pol);
        fprintf(outfile,"%.2f\t; Fraction of Situations\n",LR.conf);
        fprintf(outfile,"%.2f\t; Fraction of Time\n",LR.rel);
        fprintf(outfile,"%.2f\t; Transmitter Effective Radiated Power in Watts or dBm (optional)\n",LR.erp);
        fprintf(outfile,"\nPlease consult SPLAT! documentation for the meaning and use of this data.\n");

        fclose(outfile);

        return_value=1;

        fprintf(stderr,"\n\n%c*** There were problems reading your \"%s\" file! ***\nA \"splat.lrp\" file was written to your directory with default data.\n",7,filename);
    }

    else if (forced_read==0)
        return_value=0;

    if (forced_read && (fd==NULL || ok==0))
    {
        LR.eps_dielect=15.0;
        LR.sgm_conductivity=0.005;
        LR.eno_ns_surfref=301.0;
        LR.frq_mhz=300.0;
        LR.radio_climate=5;
        LR.pol=0;
        LR.conf=0.50;
        LR.rel=0.50;
        LR.erp=0.0;

        fprintf(stderr,"Default parameters have been assumed for this analysis.\n");

        return_value=1;
    }

    free(filename);

    return return_value;
}

/* Analyzes the path between the source and destination locations and
 * determines which points along the path have line-of-sight visibility
 * to the source. Points along with path having line-of-sight visibility
 * to the source at an AGL altitude equal to that of the destination
 * location are stored by setting bit 1 in the mask[][] array, which
 * are displayed in green when PPM maps are later generated by SPLAT!.
 */
int PlotPath(Site source, Site destination, char mask_value)
{
    char block;
    int x, y;
    double cos_xmtr_angle, cos_test_angle, test_alt;
    double distance, rx_alt, tx_alt;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return -1;
    }
    ReadPath(source,destination,path);

    for (y=0; y<path->length; y++)
    {
        /* Test this point only if it hasn't been already
           tested and found to be free of obstructions. */

        if ((GetMask(path->lat[y],path->lon[y])&mask_value)==0)
        {
            distance=5280.0*path->distance[y];

            if (source.amsl_flag)
                tx_alt=earthradius+source.alt;
            else
                tx_alt=earthradius+source.alt+path->elevation[0];

            /***
              if (destination.amsl_flag)
              rx_alt=earthradius+destination.alt;
              else
             ***/
            rx_alt=earthradius+destination.alt+path->elevation[y];

            /* Calculate the cosine of the elevation of the
               transmitter as seen at the temp rx point. */

            cos_xmtr_angle=((rx_alt*rx_alt)+(distance*distance)-(tx_alt*tx_alt))/(2.0*rx_alt*distance);

            for (x=y, block=0; x>=0 && block==0; x--)
            {
                distance=5280.0*(path->distance[y]-path->distance[x]);
                test_alt=earthradius+(path->elevation[x]==0.0?path->elevation[x]:path->elevation[x]+clutter);

                cos_test_angle=((rx_alt*rx_alt)+(distance*distance)-(test_alt*test_alt))/(2.0*rx_alt*distance);

                /* Compare these two angles to determine if
                   an obstruction exists.  Since we're comparing
                   the cosines of these angles rather than
                   the angles themselves, the following "if"
                   statement is reversed from what it would
                   be if the actual angles were compared. */

                if (cos_xmtr_angle>=cos_test_angle)
                    block=1;
            }

            if (block==0)
                OrMask(path->lat[y],path->lon[y],mask_value);
        }
    }

    DestroyPath(path);
    return 0;
}

/* Plots the RF path loss between source and destination points based on the
 * ITM/ITWOM propagation model, taking into account antenna pattern data if
 * available.
 */
int PlotLRPath(Site source, Site destination, unsigned char mask_value, FILE *fd)
{
    int	x, y, ifs, ofs, errnum;
    char	block=0, strmode[100];
    double	loss, azimuth, pattern=0.0, xmtr_alt,
            dest_alt, xmtr_alt2, dest_alt2, cos_rcvr_angle,
            cos_test_angle=0.0, test_alt, elevation=0.0,
            distance=0.0, four_thirds_earth, rxp, dBm,
            field_strength=0.0;

    Site temp;
    char textout[MAX_LINE_LEN];
    size_t textlen = 0;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return -1;
    }
    ReadPath(source,destination,path);

    /* XXX debug */
    totalpaths++;
    UPDATE_RUNNING_AVG(avgpathlen, path->length, totalpaths);

    /* XXX why +10? should it just be +2? Better yet, path->length+2? */
    elev_t *elev = (elev_t*)calloc(path->arysize + 10, sizeof(elev_t));

    four_thirds_earth=FOUR_THIRDS*EARTHRADIUS;

    /* Copy elevations plus clutter along path into the elev[] array. */

    for (x=1; x<path->length-1; x++)
        elev[x+2]=(path->elevation[x]==0.0?(elev_t)(path->elevation[x]*METERS_PER_FOOT):(elev_t)((clutter+path->elevation[x])*METERS_PER_FOOT));

    /* Copy ending points without clutter */

    elev[2]=(elev_t)(path->elevation[0]*METERS_PER_FOOT);
    elev[path->length+1]=(elev_t)(path->elevation[path->length-1]*METERS_PER_FOOT);

    /* Since the only energy the propagation model considers
       reaching the destination is based on what is scattered
       or deflected from the first obstruction along the path,
       we first need to find the location and elevation angle
       of that first obstruction (if it exists).  This is done
       using a 4/3rds Earth radius to match the radius used by
       the irregular terrain propagation model.  This information
       is required for properly integrating the antenna's elevation
       pattern into the calculation for overall path loss. */

    for (y=2; (y<(path->length-1) && path->distance[y]<=max_range); y++)
    {
        /* Process this point only if it
           has not already been processed. */

        if ((GetMask(path->lat[y],path->lon[y])&248)!=(mask_value<<3))
        {
            distance=5280.0*path->distance[y];

            /***
              if (source.amsl_flag)
              xmtr_alt=four_thirds_earth+source.alt;
              else
             ***/
            xmtr_alt=four_thirds_earth+source.alt+path->elevation[0];

            /***/
            if (destination.amsl_flag)
                dest_alt=four_thirds_earth+destination.alt;
            else
                /***/
                dest_alt=four_thirds_earth+destination.alt+path->elevation[y];

            dest_alt2=dest_alt*dest_alt;
            xmtr_alt2=xmtr_alt*xmtr_alt;

            /* Calculate the cosine of the elevation of
               the receiver as seen by the transmitter. */

            cos_rcvr_angle=((xmtr_alt2)+(distance*distance)-(dest_alt2))/(2.0*xmtr_alt*distance);

            if (cos_rcvr_angle>1.0)
                cos_rcvr_angle=1.0;

            if (cos_rcvr_angle<-1.0)
                cos_rcvr_angle=-1.0;

            if (got_elevation_pattern || fd!=NULL)
            {
                /* Determine the elevation angle to the first obstruction
                   along the path IF elevation pattern data is available
                   or an output (.ano) file has been designated. */

                for (x=2, block=0; (x<y && block==0); x++)
                {
                    distance=5280.0*path->distance[x];

                    test_alt=four_thirds_earth+(path->elevation[x]==0.0?path->elevation[x]:path->elevation[x]+clutter);

                    /* Calculate the cosine of the elevation
                       angle of the terrain (test point)
                       as seen by the transmitter. */

                    cos_test_angle=((xmtr_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*xmtr_alt*distance);

                    if (cos_test_angle>1.0)
                        cos_test_angle=1.0;

                    if (cos_test_angle<-1.0)
                        cos_test_angle=-1.0;

                    /* Compare these two angles to determine if
                       an obstruction exists.  Since we're comparing
                       the cosines of these angles rather than
                       the angles themselves, the sense of the
                       following "if" statement is reversed from
                       what it would be if the angles themselves
                       were compared. */

                    if (cos_rcvr_angle>=cos_test_angle)
                        block=1;
                }

                if (block)
                    elevation=((acos(cos_test_angle))/DEG2RAD)-90.0;
                else
                    elevation=((acos(cos_rcvr_angle))/DEG2RAD)-90.0;
            }

            /* Determine attenuation for each point along
               the path using ITWOM's point_to_point mode
               starting at y=2 (number_of_points = 1), the
               shortest distance terrain can play a role in
               path loss. */

            elev[0]=(elev_t)(y-1);  /* (number of points - 1) */

            /* Distance between elevation samples */

            elev[1]=(elev_t)(METERS_PER_MILE*(path->distance[y]-path->distance[y-1]));

            if (itwom)
                point_to_point(elev,source.alt*METERS_PER_FOOT,
                               destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                               LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                               LR.radio_climate, LR.pol, LR.conf, LR.rel, &loss,
                               strmode, &errnum);
            else
                point_to_point_ITM(elev,source.alt*METERS_PER_FOOT, 
                                   destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                                   LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                                   LR.radio_climate, LR.pol, LR.conf, LR.rel, &loss,
                                   strmode, &errnum);

            temp.lat=path->lat[y];
            temp.lon=path->lon[y];

            azimuth=(Azimuth(source,temp));

            if (fd!=NULL) {
                textlen = snprintf(textout, MAX_LINE_LEN, "%.7f, %.7f, %.3f, %.3f, ",
                                   path->lat[y], path->lon[y], azimuth, elevation);
            }

            /* If ERP==0, write path loss to alphanumeric
               output file.  Otherwise, write field strength
               or received power level (below), as appropriate. */

            if (fd!=NULL && LR.erp==0.0) {
                textlen += snprintf(textout + textlen, MAX_LINE_LEN - textlen, "%.2f",loss);
            }

            /* Integrate the antenna's radiation
               pattern into the overall path loss. */

            x=(int)rint(10.0*(10.0-elevation));

            if (x>=0 && x<=1000)
            {
                azimuth=rint(azimuth);

                pattern=(double)LR.antenna_pattern[(int)azimuth][x];

                if (pattern > 0.0) {
                    pattern=20.0*log10(pattern);
                    loss-=pattern;
                } else if (pattern != NO_ANTENNA_DATA) {
                    loss-=9999;
                }
            }

            if (LR.erp!=0.0)
            {
                if (dbm)
                {
                    /* dBm is based on EIRP (ERP + 2.14) */

                    rxp=LR.erp/(pow(10.0,(loss-2.14)/10.0));

                    dBm=10.0*(log10(rxp*1000.0));

                    if (fd!=NULL) {
                        textlen += snprintf(textout + textlen, MAX_LINE_LEN - textlen, "%.3f",dBm);
                    }

                    /* Scale roughly between 0 and 255 */

                    ifs=200+(int)rint(dBm);

                    if (ifs<0)
                        ifs=0;

                    if (ifs>255)
                        ifs=255;

                    ofs=GetSignal(path->lat[y],path->lon[y]);

                    if (ofs>ifs)
                        ifs=ofs;

                    // writes to dem
                    PutSignal(path->lat[y],path->lon[y],(unsigned char)ifs);
                }

                else
                {
                    field_strength=(139.4+(20.0*log10(LR.frq_mhz))-loss)+(10.0*log10(LR.erp/1000.0));

                    ifs=100+(int)rint(field_strength);

                    if (ifs<0)
                        ifs=0;

                    if (ifs>255)
                        ifs=255;

                    ofs=GetSignal(path->lat[y],path->lon[y]);

                    if (ofs>ifs)
                        ifs=ofs;

                    PutSignal(path->lat[y],path->lon[y],(unsigned char)ifs);

                    if (fd!=NULL) {
                        textlen += snprintf(textout + textlen, MAX_LINE_LEN - textlen, "%.3f",field_strength);
                    }
                }
            }

            else
            {
                if (loss>255)
                    ifs=255;
                else
                    ifs=(int)rint(loss);

                ofs=GetSignal(path->lat[y],path->lon[y]);

                if (ofs<ifs && ofs!=0)
                    ifs=ofs;

                PutSignal(path->lat[y],path->lon[y],(unsigned char)ifs);
            }

            if (fd!=NULL) {
                textlen += snprintf(textout + textlen, MAX_LINE_LEN - textlen, "%s",
                                    block ? " *\n" : "\n");
            }

            fwrite(textout, 1, textlen, fd);

            /* Mark this point as having been analyzed */

            PutMask(path->lat[y],path->lon[y],(GetMask(path->lat[y],path->lon[y])&7)+(mask_value<<3));
        }
    }

    free(elev);
    DestroyPath(path);

    return 0;
}


/* Performs a 360 degree sweep around the transmitter site (source location), and
 * plots the line-of-sight coverage of the transmitter on the SPLAT! generated
 * topographic map based on a receiver located at the specified altitude (in feet
 * AGL). Results are stored in memory, and written out in the form of a topographic
 * map when the WriteImage() function is later invoked.
 */
void PlotLOSMap(Site source, double altitude, bool multithread)
{
    int y, z, count;
    Site edge;
    unsigned char x;
    double lat, lon, minwest, maxnorth, th;
    static unsigned char mask_value=1;
    char symbol[4] = {'.', 'o', 'O', 'o' };

    minwest=dpp+(double)min_west;
    maxnorth=(double)max_north-dpp;

    count=0;

    fprintf(stdout,"\nComputing line-of-sight coverage of \"%s\" with an RX antenna\nat %.2f %s AGL",
            source.name, (metric?altitude*METERS_PER_FOOT:altitude), (metric?"meters":"feet"));

    if (clutter>0.0)
        fprintf(stdout," and %.2f %s of ground clutter",
                (metric?clutter*METERS_PER_FOOT:clutter), (metric?"meters":"feet"));


    /* th=pixels/degree divided by 64 loops per
       progress indicator symbol (.oOo) printed. */

    th=ppd/64.0;

    z=(int)(th*ReduceAngle(max_west-min_west));

    fprintf(stdout, "\n\n");

    WorkQueue wq;

    if (verbose) {
        if (multithread) {
            fprintf(stdout,"Using %d threads...\n\n", wq.maxWorkers());
        }
        fprintf(stdout," 0%c to  25%c ",37,37);
        fflush(stdout);
    }

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=max_north;
        edge.lon=lon;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotPath, source, edge, mask_value));
        } else {
            PlotPath(source,edge,mask_value);
        }

        if (verbose) {
            if (++count==z) 
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    if (verbose) {
        count=0;
        fprintf(stdout,"\n25%c to  50%c ",37,37);
        fflush(stdout);
    }

    z=(int)(th*(double)(max_north-min_north));

    for (lat=maxnorth, x=0, y=0; lat>=(double)min_north; y++, lat=maxnorth-(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=min_west;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotPath, source, edge, mask_value));
        } else {
            PlotPath(source,edge,mask_value);
        }

        if (verbose) {
            if (++count==z) 
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    if (verbose) {
        count=0;
        fprintf(stdout,"\n50%c to  75%c ",37,37);
        fflush(stdout);
    }

    z=(int)(th*ReduceAngle(max_west-min_west));

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=min_north;
        edge.lon=lon;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotPath, source, edge, mask_value));
        } else {
            PlotPath(source,edge,mask_value);
        }

        if (verbose) {
            if (++count==z)
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    if (verbose) {
        count=0;
        fprintf(stdout,"\n75%c to 100%c ",37,37);
        fflush(stdout);
    }

    z=(int)(th*(double)(max_north-min_north));

    for (lat=(double)min_north, x=0, y=0; lat<(double)max_north; y++, lat=(double)min_north+(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=max_west;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotPath, source, edge, mask_value));
        } else {
            PlotPath(source,edge,mask_value);
        }

        if (verbose) {
            if (++count==z)
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    wq.waitForCompletion();

    if (verbose) {
        fprintf(stdout,"\nDone!\n");
        fflush(stdout);
    }

    /* Assign next mask value */

    switch (mask_value)
    {
        case 1:
            mask_value=8;
            break;

        case 8:
            mask_value=16;
            break;

        case 16:
            mask_value=32;
    }
}

/* Performs a 360 degree sweep around the transmitter site (source location),
 * and plots the ITM/ITWOM attenuation on the SPLAT! generated topographic map
 * based on a receiver located at the specified altitude (in feet AGL).  Results
 * are stored in memory, and written out in the form of a topographic map when
 * the WriteImageLR() or WriteImageSS() functions are later invoked.
 */
void PlotLRMap(Site source, double altitude, char *plo_filename, bool multithread)
{
    int y, z, count;
    Site edge;
    double lat, lon, minwest, maxnorth, th;
    unsigned char x;
    static unsigned char mask_value=1;
    FILE *fd=NULL;
    char symbol[4] = {'.', 'o', 'O', 'o' };

    minwest=dpp+(double)min_west;
    maxnorth=(double)max_north-dpp;

    count=0;

    if (itwom)
        fprintf(stdout,"\nComputing ITWOM ");
    else
        fprintf(stdout,"\nComputing ITM ");

    if (LR.erp==0.0)
        fprintf(stdout,"path loss");
    else
    {
        if (dbm)
            fprintf(stdout,"signal power level");
        else
            fprintf(stdout,"field strength");
    }

    fprintf(stdout," contours of \"%s\"\nout to a radius of %.2f %s with an RX antenna at %.2f %s AGL",
            source.name, (metric?max_range*KM_PER_MILE:max_range), (metric?"kilometers":"miles"),
            (metric?altitude*METERS_PER_FOOT:altitude),(metric?"meters":"feet"));

    if (clutter>0.0)
        fprintf(stdout,"\nand %.2f %s of ground clutter",
                (metric?clutter*METERS_PER_FOOT:clutter), (metric?"meters":"feet"));

    if (plo_filename[0]!=0)
        fd=fopen(plo_filename,"wb");

    if (fd!=NULL)
    {
        /* Write header information to output file */

        fprintf(fd,"%d, %d\t; max_west, min_west\n%d, %d\t; max_north, min_north\n",max_west, min_west, max_north, min_north);
    }

    /* th=pixels/degree divided by 64 loops per
       progress indicator symbol (.oOo) printed. */

    th=ppd/64.0;

    z=(int)(th*ReduceAngle(max_west-min_west));

    fprintf(stdout, "\n\n");

    WorkQueue wq;

    if (verbose) {
        if (multithread) {
            fprintf(stdout,"Using %d threads...\n\n", wq.maxWorkers());
        }
        fprintf(stdout," 0%c to  25%c ",37,37);
        fflush(stdout);
    }

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=max_north;
        edge.lon=lon;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotLRPath, source, edge, mask_value, fd));
        } else {
            PlotLRPath(source,edge,mask_value,fd);
        }

        if (verbose) {
            if (++count==z) 
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    if (verbose) {
        count=0;
        fprintf(stdout,"\n25%c to  50%c ",37,37);
        fflush(stdout);
    }

    z=(int)(th*(double)(max_north-min_north));

    for (lat=maxnorth, x=0, y=0; lat>=(double)min_north; y++, lat=maxnorth-(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=min_west;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotLRPath, source, edge, mask_value, fd));
        } else {
            PlotLRPath(source,edge,mask_value,fd);
        }

        if (verbose) {
            if (++count==z) 
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    if (verbose) {
        count=0;
        fprintf(stdout,"\n50%c to  75%c ",37,37);
        fflush(stdout);
    }

    z=(int)(th*ReduceAngle(max_west-min_west));

    for (lon=minwest, x=0, y=0; (LonDiff(lon,(double)max_west)<=0.0); y++, lon=minwest+(dpp*(double)y))
    {
        if (lon>=360.0)
            lon-=360.0;

        edge.lat=min_north;
        edge.lon=lon;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotLRPath, source, edge, mask_value, fd));
        } else {
            PlotLRPath(source,edge,mask_value,fd);
        }

        if (verbose) {
            if (++count==z)
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    if (verbose) {
        count=0;
        fprintf(stdout,"\n75%c to 100%c ",37,37);
        fflush(stdout);
    }

    z=(int)(th*(double)(max_north-min_north));

    for (lat=(double)min_north, x=0, y=0; lat<(double)max_north; y++, lat=(double)min_north+(dpp*(double)y))
    {
        edge.lat=lat;
        edge.lon=max_west;
        edge.alt=altitude;

        if (multithread) {
            wq.submit(std::bind(PlotLRPath, source, edge, mask_value, fd));
        } else {
            PlotLRPath(source,edge,mask_value,fd);
        }

        if (verbose) {
            if (++count==z)
            {
                fprintf(stdout,"%c",symbol[x]);
                fflush(stdout);
                count=0;

                if (x==3)
                    x=0;
                else
                    x++;
            }
        }
    }

    wq.waitForCompletion();

    if (fd!=NULL)
        fclose(fd);

    if (verbose) {
        fprintf(stdout,"\nDone!\n");
        fprintf(stdout,"\nThere were %d paths with an average length of %d elements.\n", totalpaths, (int)avgpathlen);
        fflush(stdout);
    }

    if (mask_value<30)
        mask_value++;
}

/* Load signal-strength color values for use when plotting maps.
 */
void LoadSignalColors(Site xmtr)
{
    int x, y, ok, val[4];
    char *filename, buf[MAX_LINE_LEN], *pointer=NULL;
    FILE *fd=NULL;

    filename = copyFilename(xmtr.filename, "scf");

    /* Default values */

    region.level[0]=128;
    region.color[0][0]=255;
    region.color[0][1]=0;
    region.color[0][2]=0;

    region.level[1]=118;
    region.color[1][0]=255;
    region.color[1][1]=165;
    region.color[1][2]=0;

    region.level[2]=108;
    region.color[2][0]=255;
    region.color[2][1]=206;
    region.color[2][2]=0;

    region.level[3]=98;
    region.color[3][0]=255;
    region.color[3][1]=255;
    region.color[3][2]=0;

    region.level[4]=88;
    region.color[4][0]=184;
    region.color[4][1]=255;
    region.color[4][2]=0;

    region.level[5]=78;
    region.color[5][0]=0;
    region.color[5][1]=255;
    region.color[5][2]=0;

    region.level[6]=68;
    region.color[6][0]=0;
    region.color[6][1]=208;
    region.color[6][2]=0;

    region.level[7]=58;
    region.color[7][0]=0;
    region.color[7][1]=196;
    region.color[7][2]=196;

    region.level[8]=48;
    region.color[8][0]=0;
    region.color[8][1]=148;
    region.color[8][2]=255;

    region.level[9]=38;
    region.color[9][0]=80;
    region.color[9][1]=80;
    region.color[9][2]=255;

    region.level[10]=28;
    region.color[10][0]=0;
    region.color[10][1]=38;
    region.color[10][2]=255;

    region.level[11]=18;
    region.color[11][0]=142;
    region.color[11][1]=63;
    region.color[11][2]=255;

    region.level[12]=8;
    region.color[12][0]=140;
    region.color[12][1]=0;
    region.color[12][2]=128;

    region.levels=13;

    fd=fopen("splat.scf","r");

    if (fd==NULL)
        fd=fopen(filename,"r");

    if (fd==NULL)
    {
        fd=fopen(filename,"w");

        fprintf(fd,"; SPLAT! Auto-generated Signal Color Definition (\"%s\") File\n",filename);
        fprintf(fd,";\n; Format for the parameters held in this file is as follows:\n;\n");
        fprintf(fd,";	dBuV/m: red, green, blue\n;\n");
        fprintf(fd,"; ...where \"dBuV/m\" is the signal strength (in dBuV/m) and\n");
        fprintf(fd,"; \"red\", \"green\", and \"blue\" are the corresponding RGB color\n");
        fprintf(fd,"; definitions ranging from 0 to 255 for the region specified.\n");
        fprintf(fd,";\n; The following parameters may be edited and/or expanded\n");
        fprintf(fd,"; for future runs of SPLAT!  A total of 32 contour regions\n");
        fprintf(fd,"; may be defined in this file.\n;\n;\n");

        for (x=0; x<region.levels; x++)
            fprintf(fd,"%3d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

        fclose(fd);
    }
    else
    {
        x=0;
        fgets(buf,MAX_LINE_LEN,fd);

        while (x<32 && feof(fd)==0)
        {
            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

            if (ok==4)
            {
                for (y=0; y<4; y++)
                {
                    if (val[y]>255)
                        val[y]=255;

                    if (val[y]<0)
                        val[y]=0;
                }

                region.level[x]=val[0];
                region.color[x][0]=val[1];
                region.color[x][1]=val[2];
                region.color[x][2]=val[3];
                x++;
            }

            fgets(buf,MAX_LINE_LEN,fd);
        }

        fclose(fd);
        region.levels=x;
    }

    free(filename);
}

/* Load signal-loss color values for use when plotting maps.
 */
void LoadLossColors(Site xmtr)
{
    int x, y, ok, val[4];
    char *filename, buf[MAX_LINE_LEN], *pointer=NULL;
    FILE *fd=NULL;

    filename = copyFilename(xmtr.filename, "lcf");

    /* Default values */

    region.level[0]=80;
    region.color[0][0]=255;
    region.color[0][1]=0;
    region.color[0][2]=0;

    region.level[1]=90;
    region.color[1][0]=255;
    region.color[1][1]=128;
    region.color[1][2]=0;

    region.level[2]=100;
    region.color[2][0]=255;
    region.color[2][1]=165;
    region.color[2][2]=0;

    region.level[3]=110;
    region.color[3][0]=255;
    region.color[3][1]=206;
    region.color[3][2]=0;

    region.level[4]=120;
    region.color[4][0]=255;
    region.color[4][1]=255;
    region.color[4][2]=0;

    region.level[5]=130;
    region.color[5][0]=184;
    region.color[5][1]=255;
    region.color[5][2]=0;

    region.level[6]=140;
    region.color[6][0]=0;
    region.color[6][1]=255;
    region.color[6][2]=0;

    region.level[7]=150;
    region.color[7][0]=0;
    region.color[7][1]=208;
    region.color[7][2]=0;

    region.level[8]=160;
    region.color[8][0]=0;
    region.color[8][1]=196;
    region.color[8][2]=196;

    region.level[9]=170;
    region.color[9][0]=0;
    region.color[9][1]=148;
    region.color[9][2]=255;

    region.level[10]=180;
    region.color[10][0]=80;
    region.color[10][1]=80;
    region.color[10][2]=255;

    region.level[11]=190;
    region.color[11][0]=0;
    region.color[11][1]=38;
    region.color[11][2]=255;

    region.level[12]=200;
    region.color[12][0]=142;
    region.color[12][1]=63;
    region.color[12][2]=255;

    region.level[13]=210;
    region.color[13][0]=196;
    region.color[13][1]=54;
    region.color[13][2]=255;

    region.level[14]=220;
    region.color[14][0]=255;
    region.color[14][1]=0;
    region.color[14][2]=255;

    region.level[15]=230;
    region.color[15][0]=255;
    region.color[15][1]=194;
    region.color[15][2]=204;

    region.levels=16;

    fd=fopen("splat.lcf","r");

    if (fd==NULL)
        fd=fopen(filename,"r");

    if (fd==NULL)
    {
        fd=fopen(filename,"w");

        fprintf(fd,"; SPLAT! Auto-generated Path-Loss Color Definition (\"%s\") File\n",filename);
        fprintf(fd,";\n; Format for the parameters held in this file is as follows:\n;\n");
        fprintf(fd,";	dB: red, green, blue\n;\n");
        fprintf(fd,"; ...where \"dB\" is the path loss (in dB) and\n");
        fprintf(fd,"; \"red\", \"green\", and \"blue\" are the corresponding RGB color\n");
        fprintf(fd,"; definitions ranging from 0 to 255 for the region specified.\n");
        fprintf(fd,";\n; The following parameters may be edited and/or expanded\n");
        fprintf(fd,"; for future runs of SPLAT!  A total of 32 contour regions\n");
        fprintf(fd,"; may be defined in this file.\n;\n;\n");

        for (x=0; x<region.levels; x++)
            fprintf(fd,"%3d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

        fclose(fd);
    }

    else
    {
        x=0;
        fgets(buf,MAX_LINE_LEN,fd);

        while (x<32 && feof(fd)==0)
        {
            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

            if (ok==4)
            {
                for (y=0; y<4; y++)
                {
                    if (val[y]>255)
                        val[y]=255;

                    if (val[y]<0)
                        val[y]=0;
                }

                region.level[x]=val[0];
                region.color[x][0]=val[1];
                region.color[x][1]=val[2];
                region.color[x][2]=val[3];
                x++;
            }

            fgets(buf,MAX_LINE_LEN,fd);
        }

        fclose(fd);
        region.levels=x;
    }

    free(filename);
}

/* Load signal-strength (dBm) color values for use when plotting maps.
 */
void LoadDBMColors(Site xmtr)
{
    int x, y, ok, val[4];
    char *filename, buf[MAX_LINE_LEN], *pointer=NULL;
    FILE *fd=NULL;

    filename = copyFilename(xmtr.filename, "dcf");

    /* Default values */

    region.level[0]=0;
    region.color[0][0]=255;
    region.color[0][1]=0;
    region.color[0][2]=0;

    region.level[1]=-10;
    region.color[1][0]=255;
    region.color[1][1]=128;
    region.color[1][2]=0;

    region.level[2]=-20;
    region.color[2][0]=255;
    region.color[2][1]=165;
    region.color[2][2]=0;

    region.level[3]=-30;
    region.color[3][0]=255;
    region.color[3][1]=206;
    region.color[3][2]=0;

    region.level[4]=-40;
    region.color[4][0]=255;
    region.color[4][1]=255;
    region.color[4][2]=0;

    region.level[5]=-50;
    region.color[5][0]=184;
    region.color[5][1]=255;
    region.color[5][2]=0;

    region.level[6]=-60;
    region.color[6][0]=0;
    region.color[6][1]=255;
    region.color[6][2]=0;

    region.level[7]=-70;
    region.color[7][0]=0;
    region.color[7][1]=208;
    region.color[7][2]=0;

    region.level[8]=-80;
    region.color[8][0]=0;
    region.color[8][1]=196;
    region.color[8][2]=196;

    region.level[9]=-90;
    region.color[9][0]=0;
    region.color[9][1]=148;
    region.color[9][2]=255;

    region.level[10]=-100;
    region.color[10][0]=80;
    region.color[10][1]=80;
    region.color[10][2]=255;

    region.level[11]=-110;
    region.color[11][0]=0;
    region.color[11][1]=38;
    region.color[11][2]=255;

    region.level[12]=-120;
    region.color[12][0]=142;
    region.color[12][1]=63;
    region.color[12][2]=255;

    region.level[13]=-130;
    region.color[13][0]=196;
    region.color[13][1]=54;
    region.color[13][2]=255;

    region.level[14]=-140;
    region.color[14][0]=255;
    region.color[14][1]=0;
    region.color[14][2]=255;

    region.level[15]=-150;
    region.color[15][0]=255;
    region.color[15][1]=194;
    region.color[15][2]=204;

    region.levels=16;

    fd=fopen("splat.dcf","r");

    if (fd==NULL)
        fd=fopen(filename,"r");

    if (fd==NULL)
    {
        fd=fopen(filename,"w");

        fprintf(fd,"; SPLAT! Auto-generated DBM Signal Level Color Definition (\"%s\") File\n",filename);
        fprintf(fd,";\n; Format for the parameters held in this file is as follows:\n;\n");
        fprintf(fd,";	dBm: red, green, blue\n;\n");
        fprintf(fd,"; ...where \"dBm\" is the received signal power level between +40 dBm\n");
        fprintf(fd,"; and -200 dBm, and \"red\", \"green\", and \"blue\" are the corresponding\n");
        fprintf(fd,"; RGB color definitions ranging from 0 to 255 for the region specified.\n");
        fprintf(fd,";\n; The following parameters may be edited and/or expanded\n");
        fprintf(fd,"; for future runs of SPLAT!  A total of 32 contour regions\n");
        fprintf(fd,"; may be defined in this file.\n;\n;\n");

        for (x=0; x<region.levels; x++)
            fprintf(fd,"%+4d: %3d, %3d, %3d\n",region.level[x], region.color[x][0], region.color[x][1], region.color[x][2]);

        fclose(fd);
    }

    else
    {
        x=0;
        fgets(buf,MAX_LINE_LEN,fd);

        while (x<32 && feof(fd)==0)
        {
            pointer=strchr(buf,';');

            if (pointer!=NULL)
                *pointer=0;

            ok=sscanf(buf,"%d: %d, %d, %d", &val[0], &val[1], &val[2], &val[3]);

            if (ok==4)
            {
                if (val[0]<-200)
                    val[0]=-200;

                if (val[0]>+40)
                    val[0]=+40;

                region.level[x]=val[0];

                for (y=1; y<4; y++)
                {
                    if (val[y]>255)
                        val[y]=255;

                    if (val[y]<0)
                        val[y]=0;
                }

                region.color[x][0]=val[1];
                region.color[x][1]=val[2];
                region.color[x][2]=val[3];
                x++;
            }

            fgets(buf,MAX_LINE_LEN,fd);
        }

        fclose(fd);
        region.levels=x;
    }

    free(filename);
}

/* Generates a topographic map image based on logarithmically scaled
 * topology data as well as the content of flags held in the mask[][] array.
 * The image created is rotated counter-clockwise 90 degrees from its
 * representation in DEM[][] so that north points up and east points right.
 */
void WriteImage(char *filename, ImageType imagetype, bool geo, bool kml, bool ngs, Site *xmtr, unsigned int txsites)
{
    const char *suffix;
    char *mapfile, *geofile, *kmlfile;
#if DO_KMZ
    char *kmzfile;
#endif
    unsigned char mask;
    unsigned width, height, terrain;
    int x0=0, y0=0;
    DEM *dem;
    double conversion, one_over_gamma,
           north, south, east, west, minwest;
    FILE *fd;

    Pixel pixel = 0;
    ImageWriter iw;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    if (!filename || filename[0]=='\0')
    {
        filename = xmtr[0].filename;
    }

    switch (imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            suffix = "png";
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            suffix = "jpg";
            break;
#endif
        case IMAGETYPE_PPM:
            suffix = "ppm";
            break;
    }
    mapfile = copyFilename(filename, suffix);
    geofile = copyFilename(filename, "geo");
    kmlfile = copyFilename(filename, "kml");
#if DO_KMZ
    kmzfile = copyFilename(filename, "kmz");
#endif

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;
    south=(double)min_north;
    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);

    if (!kml && geo)
    {
        fd=fopen(geofile,"wb");

        fprintf(fd,"FILENAME\t%s\n",mapfile);
        fprintf(fd,"#\t\tX\tY\tLong\t\tLat\n");
        fprintf(fd,"TIEPOINT\t0\t0\t%.3f\t\t%.3f\n",west,north);
        fprintf(fd,"TIEPOINT\t%u\t%u\t%.3f\t\t%.3f\n",width-1,height-1,east,south);
        fprintf(fd,"IMAGESIZE\t%u\t%u\n",width,height);
        fprintf(fd,"#\n# Auto Generated by %s v%s\n#\n",SPLAT_NAME,SPLAT_VERSION);

        fclose(fd);
    }

    if (kml && !geo)
    {
        fd=fopen(kmlfile,"wb");

        fprintf(fd,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fd,"<kml xmlns=\"http://earth.google.com/kml/2.1\">\n");
        fprintf(fd,"  <Folder>\n");
        fprintf(fd,"   <name>%s</name>\n",SPLAT_NAME);
        fprintf(fd,"	 <description>Line-of-Sight Contour</description>\n");
        fprintf(fd,"	   <GroundOverlay>\n");
        fprintf(fd,"		  <name>%s Line-of-Sight Contour</name>\n",SPLAT_NAME);
        fprintf(fd,"		  <description>SPLAT! Coverage</description>\n");
        fprintf(fd,"		  <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		  <Icon>\n");
        fprintf(fd,"			  <href>%s</href>\n",mapfile);
        fprintf(fd,"		</Icon>\n");
        /* fprintf(fd,"			<opacity>128</opacity>\n"); */
        fprintf(fd,"			<LatLonBox>\n");
        fprintf(fd,"			   <north>%.5f</north>\n",north);
        fprintf(fd,"			   <south>%.5f</south>\n",south);
        fprintf(fd,"			   <east>%.5f</east>\n",east);
        fprintf(fd,"			   <west>%.5f</west>\n",west);
        fprintf(fd,"			   <rotation>0.0</rotation>\n");
        fprintf(fd,"			</LatLonBox>\n");
        fprintf(fd,"	   </GroundOverlay>\n");

        for (unsigned int x=0; x<txsites; x++)
        {
            fprintf(fd,"	 <Placemark>\n");
            fprintf(fd,"	   <name>%s</name>\n",xmtr[x].name);
            fprintf(fd,"	   <visibility>1</visibility>\n");
            fprintf(fd,"	   <Style>\n");
            fprintf(fd,"	   <IconStyle>\n");
            fprintf(fd,"		<Icon>\n");
            fprintf(fd,"		  <href>root://icons/palette-5.png</href>\n");
            fprintf(fd,"		  <x>224</x>\n");
            fprintf(fd,"		  <y>224</y>\n");
            fprintf(fd,"		  <w>32</w>\n");
            fprintf(fd,"		  <h>32</h>\n");
            fprintf(fd,"		</Icon>\n");
            fprintf(fd,"	   </IconStyle>\n");
            fprintf(fd,"	   </Style>\n");
            fprintf(fd,"	  <Point>\n");
            fprintf(fd,"		<extrude>1</extrude>\n");
            fprintf(fd,"		<altitudeMode>relativeToGround</altitudeMode>\n");
            fprintf(fd,"		<coordinates>%f,%f,%f</coordinates>\n",(xmtr[x].lon<180.0?-xmtr[x].lon:360.0-xmtr[x].lon), xmtr[x].lat, xmtr[x].alt);
            fprintf(fd,"	  </Point>\n");
            fprintf(fd,"	 </Placemark>\n");
        }

        fprintf(fd,"  </Folder>\n");
        fprintf(fd,"</kml>\n");

        fclose(fd);
    }

    fprintf(stdout,"Writing \"%s\" (%ux%u image)...\n",mapfile,width,height);
    fflush(stdout);

    if (ImageWriterInit(&iw,mapfile,imagetype,width,height)<0) {
        fprintf(stdout,"\nError writing \"%s\"!\n",mapfile);
        fflush(stdout);
        return;
    }

    for (int y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (int x=0, lon=max_west; x<(int)width; x++, lon=(double)max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            dem = FindDEM(lat, lon, x0, y0);
            if (dem)
            {
                mask=dem->mask[x0][y0];

                if (mask&2)
                    /* Text Labels: Red */
                    pixel=COLOR_RED;

                else if (mask&4)
                    /* County Boundaries: Light Cyan */
                    pixel=COLOR_LIGHTCYAN;

                else switch (mask&57)
                {
                    case 1:
                        /* TX1: Green */
                        pixel=COLOR_GREEN;
                        break;

                    case 8:
                        /* TX2: Cyan */
                        pixel=COLOR_CYAN;
                        break;

                    case 9:
                        /* TX1 + TX2: Yellow */
                        pixel=COLOR_YELLOW;
                        break;

                    case 16:
                        /* TX3: Medium Violet */
                        pixel=COLOR_MEDIUMVIOLET;
                        break;

                    case 17:
                        /* TX1 + TX3: Pink */
                        pixel=COLOR_PINK;
                        break;

                    case 24:
                        /* TX2 + TX3: Orange */
                        pixel=COLOR_ORANGE;
                        break;

                    case 25:
                        /* TX1 + TX2 + TX3: Dark Green */
                        pixel=COLOR_DARKGREEN;
                        break;

                    case 32:
                        /* TX4: Sienna 1 */
                        pixel=COLOR_SIENNA;
                        break;

                    case 33:
                        /* TX1 + TX4: Green Yellow */
                        pixel=COLOR_GREENYELLOW;
                        break;

                    case 40:
                        /* TX2 + TX4: Dark Sea Green 1 */
                        pixel=COLOR_DARKSEAGREEN1;
                        break;

                    case 41:
                        /* TX1 + TX2 + TX4: Blanched Almond */
                        pixel=COLOR_BLANCHEDALMOND;
                        break;

                    case 48:
                        /* TX3 + TX4: Dark Turquoise */
                        pixel=COLOR_DARKTURQUOISE;
                        break;

                    case 49:
                        /* TX1 + TX3 + TX4: Medium Spring Green */
                        pixel=COLOR_MEDIUMSPRINGGREEN;
                        break;

                    case 56:
                        /* TX2 + TX3 + TX4: Tan */
                        pixel=COLOR_TAN;
                        break;

                    case 57:
                        /* TX1 + TX2 + TX3 + TX4: Gold2 */
                        pixel=COLOR_GOLD2;
                        break;

                    default:
                        if (ngs)  /* No terrain */
                            pixel=COLOR_WHITE;
                        else
                        {
                            /* Sea-level: Medium Blue */
                            if (dem->data[x0][y0]==0)
                                pixel = COLOR_MEDIUMBLUE;
                            else
                            {
                                /* Elevation: Greyscale */
                                terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                pixel=RGB(terrain,terrain,terrain);
                            }
                        }
                }
            }
            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */

                pixel=COLOR_BLACK;
            }

            ImageWriterAppendPixel(&iw, pixel);
        }

        ImageWriterEmitLine(&iw);
    }

    ImageWriterFinish(&iw);

#if DO_KMZ
    if (kml) {
        bool success = false;
        struct zip_t *zip = zip_open(kmzfile, ZIP_DEFAULT_COMPRESSION_LEVEL, 'w');
        if (zip) {
            /* Pack the KML */
            if (zip_entry_open(zip, kmlfile) == 0) {
                if (zip_entry_fwrite(zip, kmlfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            /* Pack the map image */
            if (zip_entry_open(zip, mapfile) == 0) {
                if (zip_entry_fwrite(zip, mapfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            zip_close(zip);
        }

        if (success) {
            unlink(mapfile);
            unlink(kmlfile);
            fprintf(stdout, "\nKMZ file written to: \"%s\"\n", kmzfile);
        } else {
            unlink(kmzfile);
            fprintf(stdout, "\nCouldn't create KMZ file.\n");
        }
        free(kmzfile);
    }
#endif

    free(mapfile);
    free(geofile);
    free(kmlfile);

    fprintf(stdout,"Done!\n");
    fflush(stdout);
}

/* Generates a topographic map based on the content of flags held in the mask[][] 
 * array (only). The image created is rotated counter-clockwise 90 degrees from
 * its representation in DEM[][] so that north points up and east points right.
 */
void WriteImageLR(char *filename, ImageType imagetype, bool geo, bool kml, bool ngs, Site *xmtr, unsigned int txsites)
{
    const char *suffix;
    char *mapfile, *geofile, *kmlfile, *ckfile;
#if DO_KMZ
    char *kmzfile;
#endif
    unsigned int width, height, red, green, blue, terrain=0;
    unsigned int imgheight, imgwidth;
    unsigned char mask;
    int indx, x, y, z, colorwidth, x0, y0, loss, level,
        hundreds, tens, units, match;
    DEM *dem;
    double lat, lon, conversion, one_over_gamma,
           north, south, east, west, minwest;
    FILE *fd;

    Pixel pixel = 0;
    ImageWriter iw;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    LoadLossColors(xmtr[0]);

    if (!filename || filename[0]=='\0')
    {
        filename = xmtr[0].filename;
    }

    switch (imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            suffix = "png";
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            suffix = "jpg";
            break;
#endif
        case IMAGETYPE_PPM:
            suffix = "ppm";
            break;
    }
    mapfile = copyFilename(filename, suffix);
    geofile = copyFilename(filename, "geo");
    kmlfile = copyFilename(filename, "kml");
#if DO_KMZ
    kmzfile = copyFilename(filename, "kmz");
#endif

    /* ick */
    ckfile = copyFilename(filename, NULL);
    ckfile = (char*)realloc(ckfile, strlen(ckfile)+strlen(suffix)+5);
    strcat(ckfile, "-ck.");
    strcat(ckfile, suffix);

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    if (kml || geo)
        south=(double)min_north;	/* No bottom legend */
    else
        south=(double)min_north-(30.0/ppd); /* 30 pixels for bottom legend */

    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);

    if (geo && !kml)
    {
        fd=fopen(geofile,"wb");

        fprintf(fd,"FILENAME\t%s\n",mapfile);
        fprintf(fd,"#\t\tX\tY\tLong\t\tLat\n");
        fprintf(fd,"TIEPOINT\t0\t0\t%.3f\t\t%.3f\n",west,north);

        fprintf(fd,"TIEPOINT\t%u\t%u\t%.3f\t\t%.3f\n",width-1,height-1,east,south);
        fprintf(fd,"IMAGESIZE\t%u\t%u\n",width,height);

        fprintf(fd,"#\n# Auto Generated by %s v%s\n#\n",SPLAT_NAME,SPLAT_VERSION);

        fclose(fd);
    }

    if (kml && !geo)
    {   
        fd=fopen(kmlfile,"wb");


        fprintf(fd,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fd,"<kml xmlns=\"http://earth.google.com/kml/2.1\">\n");
        fprintf(fd,"<!-- Generated by %s Version %s -->\n",SPLAT_NAME,SPLAT_VERSION);
        fprintf(fd,"  <Folder>\n");
        fprintf(fd,"   <name>%s</name>\n",SPLAT_NAME);
        fprintf(fd,"	 <description>%s Transmitter Path Loss Overlay</description>\n",xmtr[0].name);
        fprintf(fd,"	   <GroundOverlay>\n");
        fprintf(fd,"		  <name>SPLAT! Path Loss Overlay</name>\n");
        fprintf(fd,"		  <description>SPLAT! Coverage</description>\n");
        fprintf(fd,"		  <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		  <Icon>\n");
        fprintf(fd,"			  <href>%s</href>\n",basename_s(mapfile));
        fprintf(fd,"		  </Icon>\n");
        /* fprintf(fd,"			<opacity>128</opacity>\n"); */
        fprintf(fd,"			<LatLonBox>\n");
        fprintf(fd,"			   <north>%.5f</north>\n",north);
        fprintf(fd,"			   <south>%.5f</south>\n",south);
        fprintf(fd,"			   <east>%.5f</east>\n",east);
        fprintf(fd,"			   <west>%.5f</west>\n",west);
        fprintf(fd,"			   <rotation>0.0</rotation>\n");
        fprintf(fd,"			</LatLonBox>\n");
        fprintf(fd,"	   </GroundOverlay>\n");
        fprintf(fd,"	   <ScreenOverlay>\n");
        fprintf(fd,"		  <name>Color Key</name>\n");
        fprintf(fd,"		  <description>Contour Color Key</description>\n");
        fprintf(fd,"		  <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		  <Icon>\n");
        fprintf(fd,"			<href>%s</href>\n",basename_s(ckfile));
        fprintf(fd,"		  </Icon>\n");
        fprintf(fd,"		  <overlayXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <screenXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <rotationXY x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <size x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"	   </ScreenOverlay>\n");

        for (unsigned int i=0; i<txsites; i++)
        {
            fprintf(fd,"	 <Placemark>\n");
            fprintf(fd,"	   <name>%s</name>\n",xmtr[i].name);
            fprintf(fd,"	   <visibility>1</visibility>\n");
            fprintf(fd,"	   <Style>\n");
            fprintf(fd,"	   <IconStyle>\n");
            fprintf(fd,"		<Icon>\n");
            fprintf(fd,"		  <href>root://icons/palette-5.png</href>\n");
            fprintf(fd,"		  <x>224</x>\n");
            fprintf(fd,"		  <y>224</y>\n");
            fprintf(fd,"		  <w>32</w>\n");
            fprintf(fd,"		  <h>32</h>\n");
            fprintf(fd,"		</Icon>\n");
            fprintf(fd,"	   </IconStyle>\n");
            fprintf(fd,"	   </Style>\n");
            fprintf(fd,"	  <Point>\n");
            fprintf(fd,"		<extrude>1</extrude>\n");
            fprintf(fd,"		<altitudeMode>relativeToGround</altitudeMode>\n");
            fprintf(fd,"		<coordinates>%f,%f,%f</coordinates>\n",(xmtr[i].lon<180.0?-xmtr[i].lon:360.0-xmtr[i].lon), xmtr[i].lat, xmtr[i].alt);
            fprintf(fd,"	  </Point>\n");
            fprintf(fd,"	 </Placemark>\n");
        }

        fprintf(fd,"  </Folder>\n");
        fprintf(fd,"</kml>\n");

        fclose(fd);
    }

    if (kml || geo)
    {
        /* No bottom legend */
        imgwidth = width;
        imgheight = height;
    }
    else
    {
        /* Allow space for bottom legend */
        imgwidth = width;
        imgheight = height + 30;
    }

    fprintf(stdout,"\nWriting LR map \"%s\" (%ux%u image)... ",mapfile,imgwidth,imgheight);
    fflush(stdout);

    if (ImageWriterInit(&iw,mapfile,imagetype,imgwidth,imgheight)<0) {
        fprintf(stdout,"\nError writing \"%s\"!\n",mapfile);
        fflush(stdout);
        return;
    }

    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            dem = FindDEM(lat, lon, x0, y0);
            if (dem)
            {
                mask=dem->mask[x0][y0];
                loss=(dem->signal[x0][y0]);

                match=255;

                red=0;
                green=0;
                blue=0;

                if (loss<=region.level[0])
                    match=0;
                else
                {
                    for (z=1; (z<region.levels && match==255); z++)
                    {
                        if (loss>=region.level[z-1] && loss<region.level[z])
                            match=z;
                    }
                }

                if (match<region.levels)
                {
                    if (smooth_contours && match>0)
                    {
                        red=(unsigned)interpolate(region.color[match-1][0],region.color[match][0],region.level[match-1],region.level[match],loss);
                        green=(unsigned)interpolate(region.color[match-1][1],region.color[match][1],region.level[match-1],region.level[match],loss);
                        blue=(unsigned)interpolate(region.color[match-1][2],region.color[match][2],region.level[match-1],region.level[match],loss);
                    }

                    else
                    {
                        red=region.color[match][0];
                        green=region.color[match][1];
                        blue=region.color[match][2];
                    }
                }

                if (mask&2)
                {
                    /* Text Labels: Red or otherwise */

                    if (red>=180 && green<=75 && blue<=75 && loss!=0)
                        pixel=RGB(255^red,255^green,255^blue);
                    else
                        pixel=COLOR_RED;
                }
                else if (mask&4)
                {
                    /* County Boundaries: Black */
                    pixel=COLOR_BLACK;
                }
                else
                {
                    if (loss==0 || (contour_threshold!=0 && loss>abs(contour_threshold)))
                    {
                        if (ngs)  /* No terrain */
                            pixel=COLOR_WHITE;
                        else
                        {
                            /* Display land or sea elevation */

                            if (dem->data[x0][y0]==0)
                                pixel=COLOR_MEDIUMBLUE;
                            else
                            {
                                terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                pixel=RGB(terrain,terrain,terrain);
                            }
                        }
                    }
                    else
                    {
                        /* Plot path loss in color */

                        if (red!=0 || green!=0 || blue!=0)
                            pixel=RGB(red,green,blue);

                        else  /* terrain / sea-level */
                        {
                            if (dem->data[x0][y0]==0)
                                pixel=COLOR_MEDIUMBLUE;
                            else
                            {
                                /* Elevation: Greyscale */
                                terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                pixel=RGB(terrain,terrain,terrain);
                            }
                        }
                    }
                }
            }
            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */
                pixel=COLOR_BLACK;
            }

            ImageWriterAppendPixel(&iw, pixel);
        }

        ImageWriterEmitLine(&iw);
    }

    if (!kml && !geo)
    {
        /* Display legend along bottom of image
         * if not generating .kml or .geo output.
         */

        colorwidth=(int)rint((float)width/(float)region.levels);

        for (y0=0; y0<30; y0++)
        {
            for (x0=0; x0<(int)width; x0++)
            {
                indx=x0/colorwidth;
                x=x0%colorwidth;
                level=region.level[indx];

                hundreds=level/100;

                if (hundreds>0)
                    level-=(hundreds*100);

                tens=level/10;

                if (tens>0)
                    level-=(tens*10);

                units=level;

                if (y0>=8 && y0<=23)
                {  
                    if (hundreds>0)
                    {
                        if (x>=11 && x<=18)	 
                            if (fontdata[16*(hundreds+'0')+(y0-8)]&(128>>(x-11)))
                                indx=255; 
                    }

                    if (tens>0 || hundreds>0)
                    {
                        if (x>=19 && x<=26)	 
                            if (fontdata[16*(tens+'0')+(y0-8)]&(128>>(x-19)))
                                indx=255;
                    }

                    if (x>=27 && x<=34)
                        if (fontdata[16*(units+'0')+(y0-8)]&(128>>(x-27)))
                            indx=255;

                    if (x>=42 && x<=49)
                        if (fontdata[16*('d')+(y0-8)]&(128>>(x-42)))
                            indx=255;

                    if (x>=50 && x<=57)
                        if (fontdata[16*('B')+(y0-8)]&(128>>(x-50)))
                            indx=255;
                }

                if (indx>region.levels)
                    pixel=COLOR_BLACK;
                else
                {
                    red=region.color[indx][0];
                    green=region.color[indx][1];
                    blue=region.color[indx][2];

                    pixel=RGB(red,green,blue);
                }

                ImageWriterAppendPixel(&iw, pixel);
            }

            ImageWriterEmitLine(&iw);
        }
    }

    ImageWriterFinish(&iw);

    if (kml)
    {
        /* Write colorkey image file */

        height=30*region.levels;
        width=100;

        if (ImageWriterInit(&iw,ckfile,imagetype,width,height)<0) {
            fprintf(stdout,"\nError writing \"%s\"!\n",ckfile);
            fflush(stdout);
            return;
        }

        for (y0=0; y0<(int)height; y0++)
        {
            for (x0=0; x0<(int)width; x0++)
            {
                indx=y0/30;
                x=x0;
                level=region.level[indx];

                hundreds=level/100;

                if (hundreds>0)
                    level-=(hundreds*100);

                tens=level/10;

                if (tens>0)
                    level-=(tens*10);

                units=level;

                if ((y0%30)>=8 && (y0%30)<=23)
                {  
                    if (hundreds>0)
                    {
                        if (x>=11 && x<=18)	 
                            if (fontdata[16*(hundreds+'0')+((y0%30)-8)]&(128>>(x-11)))
                                indx=255; 
                    }

                    if (tens>0 || hundreds>0)
                    {
                        if (x>=19 && x<=26)	 
                            if (fontdata[16*(tens+'0')+((y0%30)-8)]&(128>>(x-19)))
                                indx=255;
                    }

                    if (x>=27 && x<=34)
                        if (fontdata[16*(units+'0')+((y0%30)-8)]&(128>>(x-27)))
                            indx=255;

                    if (x>=42 && x<=49)
                        if (fontdata[16*('d')+((y0%30)-8)]&(128>>(x-42)))
                            indx=255;

                    if (x>=50 && x<=57)
                        if (fontdata[16*('B')+((y0%30)-8)]&(128>>(x-50)))
                            indx=255;
                }

                if (indx>region.levels) {
                    pixel=COLOR_BLACK;
                }
                else
                {
                    red=region.color[indx][0];
                    green=region.color[indx][1];
                    blue=region.color[indx][2];

                    pixel=RGB(red,green,blue);
                }

                ImageWriterAppendPixel(&iw, pixel);
            } 

            ImageWriterEmitLine(&iw);
        }

        ImageWriterFinish(&iw);

#if DO_KMZ
        bool success = false;
        struct zip_t *zip = zip_open(kmzfile, ZIP_DEFAULT_COMPRESSION_LEVEL, 'w');
        if (zip) {
            /* Pack the KML */
            if (zip_entry_open(zip, kmlfile) == 0) {
                if (zip_entry_fwrite(zip, kmlfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            /* Pack the -ck file */
            if (zip_entry_open(zip, ckfile) == 0) {
                if (zip_entry_fwrite(zip, ckfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            /* Pack the map image */
            if (zip_entry_open(zip, mapfile) == 0) {
                if (zip_entry_fwrite(zip, mapfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            zip_close(zip);
        }

        if (success) {
            unlink(mapfile);
            unlink(kmlfile);
            unlink(ckfile);
            fprintf(stdout, "\nKMZ file written to: \"%s\"\n", kmzfile);
        } else {
            unlink(kmzfile);
            fprintf(stdout, "\nCouldn't create KMZ file.\n");
        }
        free(kmzfile);
#endif
    }

    free(mapfile);
    free(geofile);
    free(kmlfile);
    free(ckfile);

    fprintf(stdout,"Done!\n");
    fflush(stdout);
}

/* Generates a topographic map based on the signal strength values held in the
 * signal[][] array. The image created is rotated counter-clockwise 90 degrees
 * from its representation in dem[][] so that north points up and east points
 * right.
 *
 * In this version of the WriteImage function the signal strength is
 *  plotted (vs the power level, plotted by WriteImageDBM).
 */
void WriteImageSS(char *filename, ImageType imagetype, bool geo, bool kml, bool ngs, Site *xmtr, unsigned int txsites)
{
    const char *suffix;
    char *mapfile, *geofile, *kmlfile, *ckfile;
#if DO_KMZ
    char *kmzfile;
#endif
    unsigned width, height, terrain, red, green, blue;
    unsigned int imgheight, imgwidth;
    unsigned char mask;
    int indx, x, y, z=1, x0, y0, signal, level, hundreds,
        tens, units, match, colorwidth;
    DEM *dem;
    double conversion, one_over_gamma, lat, lon,
           north, south, east, west, minwest;
    FILE *fd;

    Pixel pixel = 0;
    ImageWriter iw;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    LoadSignalColors(xmtr[0]);

    if (!filename || filename[0]=='\0')
    {
        filename = xmtr[0].filename;
    }

    switch (imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            suffix = "png";
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            suffix = "jpg";
            break;
#endif
        case IMAGETYPE_PPM:
            suffix = "ppm";
            break;
    }
    mapfile = copyFilename(filename, suffix);
    geofile = copyFilename(filename, "geo");
    kmlfile = copyFilename(filename, "kml");
#if DO_KMZ
    kmzfile = copyFilename(filename, "kmz");
#endif

    /* ick */
    ckfile = copyFilename(filename, NULL);
    ckfile = (char*)realloc(ckfile, strlen(ckfile)+strlen(suffix)+5);
    strcat(ckfile, "-ck.");
    strcat(ckfile, suffix);

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    if (kml || geo)
        south=(double)min_north;	/* No bottom legend */
    else
        south=(double)min_north-(30.0/ppd);	/* 30 pixels for bottom legend */

    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);

    if (geo && !kml)
    {
        fd=fopen(geofile,"wb");

        fprintf(fd,"FILENAME\t%s\n",mapfile);
        fprintf(fd,"#\t\tX\tY\tLong\t\tLat\n");
        fprintf(fd,"TIEPOINT\t0\t0\t%.3f\t\t%.3f\n",west,north);

        fprintf(fd,"TIEPOINT\t%u\t%u\t%.3f\t\t%.3f\n",width-1,height-1,east,south);
        fprintf(fd,"IMAGESIZE\t%u\t%u\n",width,height);

        fprintf(fd,"#\n# Auto Generated by %s v%s\n#\n",SPLAT_NAME,SPLAT_VERSION);

        fclose(fd);
    }

    if (kml && !geo)
    {   
        fd=fopen(kmlfile,"wb");

        fprintf(fd,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fd,"<kml xmlns=\"http://earth.google.com/kml/2.1\">\n");
        fprintf(fd,"<!-- Generated by %s Version %s -->\n",SPLAT_NAME,SPLAT_VERSION);
        fprintf(fd,"  <Folder>\n");
        fprintf(fd,"   <name>%s</name>\n",SPLAT_NAME);
        fprintf(fd,"	 <description>%s Transmitter Contours</description>\n",xmtr[0].name);
        fprintf(fd,"	   <GroundOverlay>\n");
        fprintf(fd,"		 <name>SPLAT! Signal Strength Contours</name>\n");
        fprintf(fd,"		 <description>SPLAT! Coverage</description>\n");
        fprintf(fd,"		 <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		 <Icon>\n");
        fprintf(fd,"			  <href>%s</href>\n",basename_s(mapfile));
        fprintf(fd,"		 </Icon>\n");
        /* fprintf(fd,"			<opacity>128</opacity>\n"); */
        fprintf(fd,"			<LatLonBox>\n");
        fprintf(fd,"			   <north>%.5f</north>\n",north);
        fprintf(fd,"			   <south>%.5f</south>\n",south);
        fprintf(fd,"			   <east>%.5f</east>\n",east);
        fprintf(fd,"			   <west>%.5f</west>\n",west);
        fprintf(fd,"			   <rotation>0.0</rotation>\n");
        fprintf(fd,"			</LatLonBox>\n");
        fprintf(fd,"	   </GroundOverlay>\n");
        fprintf(fd,"	   <ScreenOverlay>\n");
        fprintf(fd,"		  <name>Color Key</name>\n");
        fprintf(fd,"		  <description>Contour Color Key</description>\n");
        fprintf(fd,"		  <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		  <Icon>\n");
        fprintf(fd,"			<href>%s</href>\n",basename_s(ckfile));
        fprintf(fd,"		  </Icon>\n");
        fprintf(fd,"		  <overlayXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <screenXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <rotationXY x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <size x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"	   </ScreenOverlay>\n");

        for (unsigned int i=0; i<txsites; i++)
        {
            fprintf(fd,"	 <Placemark>\n");
            fprintf(fd,"	   <name>%s</name>\n",xmtr[i].name);
            fprintf(fd,"	   <visibility>1</visibility>\n");
            fprintf(fd,"	   <Style>\n");
            fprintf(fd,"	   <IconStyle>\n");
            fprintf(fd,"		<Icon>\n");
            fprintf(fd,"		  <href>root://icons/palette-5.png</href>\n");
            fprintf(fd,"		  <x>224</x>\n");
            fprintf(fd,"		  <y>224</y>\n");
            fprintf(fd,"		  <w>32</w>\n");
            fprintf(fd,"		  <h>32</h>\n");
            fprintf(fd,"		</Icon>\n");
            fprintf(fd,"	   </IconStyle>\n");
            fprintf(fd,"	   </Style>\n");
            fprintf(fd,"	  <Point>\n");
            fprintf(fd,"		<extrude>1</extrude>\n");
            fprintf(fd,"		<altitudeMode>relativeToGround</altitudeMode>\n");
            fprintf(fd,"		<coordinates>%f,%f,%f</coordinates>\n",(xmtr[i].lon<180.0?-xmtr[i].lon:360.0-xmtr[i].lon), xmtr[i].lat, xmtr[i].alt);
            fprintf(fd,"	  </Point>\n");
            fprintf(fd,"	 </Placemark>\n");
        }

        fprintf(fd,"  </Folder>\n");
        fprintf(fd,"</kml>\n");

        fclose(fd);
    }

    if (kml || geo)
    {
        /* No bottom legend */
        imgwidth = width;
        imgheight = height;
    }
    else
    {
        /* Allow space for bottom legend */
        imgwidth = width;
        imgheight = height + 30;
    }

    fprintf(stdout,"\nWriting Signal Strength map \"%s\" (%ux%u image)... ",mapfile,imgwidth,imgheight);
    fflush(stdout);

    if (ImageWriterInit(&iw,mapfile,imagetype,imgwidth,imgheight)<0) {
        fprintf(stdout,"\nError writing \"%s\"!\n",mapfile);
        fflush(stdout);
        return;
    }

    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            dem = FindDEM(lat, lon, x0, y0);
            if (dem)
            {
                mask=dem->mask[x0][y0];
                signal=(dem->signal[x0][y0])-100;

                match=255;

                red=0;
                green=0;
                blue=0;

                if (signal>=region.level[0])
                    match=0;
                else
                {
                    for (z=1; (z<region.levels && match==255); z++)
                    {
                        if (signal<region.level[z-1] && signal>=region.level[z])
                            match=z;
                    }
                }

                if (match<region.levels)
                {
                    if (smooth_contours && match>0)
                    {
                        red=(unsigned)interpolate(region.color[match][0],region.color[match-1][0],region.level[match],region.level[match-1],signal);
                        green=(unsigned)interpolate(region.color[match][1],region.color[match-1][1],region.level[match],region.level[match-1],signal);
                        blue=(unsigned)interpolate(region.color[match][2],region.color[match-1][2],region.level[match],region.level[match-1],signal);
                    }

                    else
                    {
                        red=region.color[match][0];
                        green=region.color[match][1];
                        blue=region.color[match][2];
                    }
                }

                if (mask&2) 
                {
                    /* Text Labels: Red or otherwise */

                    if (red>=180 && green<=75 && blue<=75)
                        pixel=RGB(255^red,255^green,255^blue);
                    else
                        pixel=COLOR_RED;
                }
                else if (mask&4)
                {
                    /* County Boundaries: Black */
                    pixel=COLOR_BLACK;
                } 
                else
                {
                    if (contour_threshold!=0 && signal<contour_threshold)
                    {
                        if (ngs)  /* No terrain */
                            pixel=COLOR_WHITE;
                        else
                        {
                            /* Display land or sea elevation */

                            if (dem->data[x0][y0]==0)
                                pixel=COLOR_MEDIUMBLUE;
                            else
                            {
                                terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                pixel=RGB(terrain,terrain,terrain);
                            }
                        }
                    }

                    else
                    {
                        /* Plot field strength regions in color */

                        if (red!=0 || green!=0 || blue!=0)
                            pixel=RGB(red,green,blue);

                        else  /* terrain / sea-level */
                        {
                            if (ngs)
                                pixel=COLOR_WHITE;
                            else
                            {
                                if (dem->data[x0][y0]==0)
                                    pixel=COLOR_MEDIUMBLUE;
                                else
                                {
                                    /* Elevation: Greyscale */
                                    terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                    pixel=RGB(terrain,terrain,terrain);
                                }
                            }
                        }
                    }
                }
            }

            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */

                pixel=COLOR_BLACK;
            }

            ImageWriterAppendPixel(&iw, pixel);
        }

        ImageWriterEmitLine(&iw);
    }

    if (!kml && !geo)
    {
        /* Display legend along bottom of image
         * if not generating .kml or .geo output.
         */

        colorwidth=(int)rint((float)width/(float)region.levels);

        for (y0=0; y0<30; y0++)
        {
            for (x0=0; x0<(int)width; x0++)
            {
                indx=x0/colorwidth;
                x=x0%colorwidth;
                level=region.level[indx];

                hundreds=level/100;

                if (hundreds>0)
                    level-=(hundreds*100);

                tens=level/10;

                if (tens>0)
                    level-=(tens*10);

                units=level;

                if (y0>=8 && y0<=23)
                {  
                    if (hundreds>0)
                    {
                        if (x>=5 && x<=12)	 
                            if (fontdata[16*(hundreds+'0')+(y0-8)]&(128>>(x-5)))
                                indx=255; 
                    }

                    if (tens>0 || hundreds>0)
                    {
                        if (x>=13 && x<=20)	 
                            if (fontdata[16*(tens+'0')+(y0-8)]&(128>>(x-13)))
                                indx=255;
                    }

                    if (x>=21 && x<=28)
                        if (fontdata[16*(units+'0')+(y0-8)]&(128>>(x-21)))
                            indx=255;

                    if (x>=36 && x<=43)
                        if (fontdata[16*('d')+(y0-8)]&(128>>(x-36)))
                            indx=255;

                    if (x>=44 && x<=51)
                        if (fontdata[16*('B')+(y0-8)]&(128>>(x-44)))
                            indx=255;

                    if (x>=52 && x<=59)
                        if (fontdata[16*(230)+(y0-8)]&(128>>(x-52)))
                            indx=255;

                    if (x>=60 && x<=67)
                        if (fontdata[16*('V')+(y0-8)]&(128>>(x-60)))
                            indx=255;

                    if (x>=68 && x<=75)
                        if (fontdata[16*('/')+(y0-8)]&(128>>(x-68)))
                            indx=255;

                    if (x>=76 && x<=83)
                        if (fontdata[16*('m')+(y0-8)]&(128>>(x-76)))
                            indx=255;
                }

                if (indx>region.levels)
                    pixel=COLOR_BLACK;
                else
                {
                    red=region.color[indx][0];
                    green=region.color[indx][1];
                    blue=region.color[indx][2];

                    pixel=RGB(red,green,blue);
                }

                ImageWriterAppendPixel(&iw, pixel);
            }

            ImageWriterEmitLine(&iw);
        }
    }

    ImageWriterFinish(&iw);

    if (kml)
    {
        /* Write colorkey image file */

        height=30*region.levels;
        width=100;

        if (ImageWriterInit(&iw,ckfile,imagetype,width,height)<0) {
            fprintf(stdout,"\nError writing \"%s\"!\n",ckfile);
            fflush(stdout);
            return;
        }

        for (y0=0; y0<(int)height; y0++)
        {
            for (x0=0; x0<(int)width; x0++)
            {
                indx=y0/30;
                x=x0;

                level=region.level[indx];

                hundreds=level/100;

                if (hundreds>0)
                    level-=(hundreds*100);

                tens=level/10;

                if (tens>0)
                    level-=(tens*10);

                units=level;

                if ((y0%30)>=8 && (y0%30)<=23)
                {  
                    if (hundreds>0)
                    {
                        if (x>=5 && x<=12)	 
                            if (fontdata[16*(hundreds+'0')+((y0%30)-8)]&(128>>(x-5)))
                                indx=255; 
                    }

                    if (tens>0 || hundreds>0)
                    {
                        if (x>=13 && x<=20)	 
                            if (fontdata[16*(tens+'0')+((y0%30)-8)]&(128>>(x-13)))
                                indx=255;
                    }

                    if (x>=21 && x<=28)
                        if (fontdata[16*(units+'0')+((y0%30)-8)]&(128>>(x-21)))
                            indx=255;

                    if (x>=36 && x<=43)
                        if (fontdata[16*('d')+((y0%30)-8)]&(128>>(x-36)))
                            indx=255;

                    if (x>=44 && x<=51)
                        if (fontdata[16*('B')+((y0%30)-8)]&(128>>(x-44)))
                            indx=255;

                    if (x>=52 && x<=59)
                        if (fontdata[16*(230)+((y0%30)-8)]&(128>>(x-52)))
                            indx=255;

                    if (x>=60 && x<=67)
                        if (fontdata[16*('V')+((y0%30)-8)]&(128>>(x-60)))
                            indx=255;

                    if (x>=68 && x<=75)
                        if (fontdata[16*('/')+((y0%30)-8)]&(128>>(x-68)))
                            indx=255;

                    if (x>=76 && x<=83)
                        if (fontdata[16*('m')+((y0%30)-8)]&(128>>(x-76)))
                            indx=255;
                }

                if (indx>region.levels) {
                    pixel=COLOR_BLACK;
                }
                else
                {
                    red=region.color[indx][0];
                    green=region.color[indx][1];
                    blue=region.color[indx][2];

                    pixel=RGB(red,green,blue);
                }

                ImageWriterAppendPixel(&iw, pixel);
            } 

            ImageWriterEmitLine(&iw);
        }

        ImageWriterFinish(&iw);

#if DO_KMZ
        bool success = false;
        struct zip_t *zip = zip_open(kmzfile, ZIP_DEFAULT_COMPRESSION_LEVEL, 'w');
        if (zip) {
            /* Pack the KML */
            if (zip_entry_open(zip, kmlfile) == 0) {
                if (zip_entry_fwrite(zip, kmlfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            /* Pack the -ck file */
            if (zip_entry_open(zip, ckfile) == 0) {
                if (zip_entry_fwrite(zip, ckfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            /* Pack the map image */
            if (zip_entry_open(zip, mapfile) == 0) {
                if (zip_entry_fwrite(zip, mapfile) == 0) {
                    success = true;
                }
                zip_entry_close(zip);
            }
            zip_close(zip);
        }

        if (success) {
            unlink(mapfile);
            unlink(kmlfile);
            unlink(ckfile);
            fprintf(stdout, "\nKMZ file written to: \"%s\"\n", kmzfile);
        } else {
            unlink(kmzfile);
            fprintf(stdout, "\nCouldn't create KMZ file.\n");
        }
        free(kmzfile);
#endif
    }

    free(mapfile);
    free(geofile);
    free(kmlfile);
    free(ckfile);

    fprintf(stdout,"Done!\n");
    fflush(stdout);
}

/* Generates a topographic map based on the signal strength values held in the
 * signal[][] array. The image created is rotated counter-clockwise 90 degrees
 * from its representation in dem[][] so that north points up and east points
 * right.
 *
 * In this version of the WriteImage function the power level is
 *  plotted (vs the signal strength, plotted by WriteImageDBM).
 */
void WriteImageDBM(char *filename, ImageType imagetype, bool geo, bool kml, bool ngs, Site *xmtr, unsigned int txsites)
{
    const char *suffix;
    char *mapfile, *geofile, *kmlfile, *ckfile;
    unsigned width, height, terrain, red, green, blue;
    unsigned int imgheight, imgwidth;
    unsigned char mask;
    DEM *dem;
    int indx, x, y, z=1, x0, y0, dBm, level, hundreds,
        tens, units, match, colorwidth;
    double conversion, one_over_gamma, lat, lon,
           north, south, east, west, minwest;
    FILE *fd;

    Pixel pixel = 0;
    ImageWriter iw;

    one_over_gamma=1.0/GAMMA;
    conversion=255.0/pow((double)(max_elevation-min_elevation),one_over_gamma);

    width=(unsigned)(ippd*ReduceAngle(max_west-min_west));
    height=(unsigned)(ippd*ReduceAngle(max_north-min_north));

    LoadDBMColors(xmtr[0]);

    if (!filename || filename[0]=='\0')
    {
        filename = xmtr[0].filename;
    }

    switch (imagetype) {
        default:
#ifdef HAVE_LIBPNG
        case IMAGETYPE_PNG:
            suffix = "png";
            break;
#endif
#ifdef HAVE_LIBJPEG
        case IMAGETYPE_JPG:
            suffix = "jpg";
            break;
#endif
        case IMAGETYPE_PPM:
            suffix = "ppm";
            break;
    }
    mapfile = copyFilename(filename, suffix);
    geofile = copyFilename(filename, "geo");
    kmlfile = copyFilename(filename, "kml");

    /* ick */
    ckfile = copyFilename(filename, NULL);
    ckfile = (char*)realloc(ckfile, strlen(ckfile)+strlen(suffix)+5);
    strcat(ckfile, "-ck.");
    strcat(ckfile, suffix);

    minwest=((double)min_west)+dpp;

    if (minwest>360.0)
        minwest-=360.0;

    north=(double)max_north-dpp;

    if (kml || geo)
        south=(double)min_north;	/* No bottom legend */
    else
        south=(double)min_north-(30.0/ppd);	/* 30 pixels for bottom legend */

    east=(minwest<180.0?-minwest:360.0-min_west);
    west=(double)(max_west<180?-max_west:360-max_west);

    if (geo && !kml)
    {
        fd=fopen(geofile,"wb");

        fprintf(fd,"FILENAME\t%s\n",mapfile);
        fprintf(fd,"#\t\tX\tY\tLong\t\tLat\n");
        fprintf(fd,"TIEPOINT\t0\t0\t%.3f\t\t%.3f\n",west,north);

        fprintf(fd,"TIEPOINT\t%u\t%u\t%.3f\t\t%.3f\n",width-1,height-1,east,south);
        fprintf(fd,"IMAGESIZE\t%u\t%u\n",width,height);

        fprintf(fd,"#\n# Auto Generated by %s v%s\n#\n",SPLAT_NAME,SPLAT_VERSION);

        fclose(fd);
    }

    if (kml && !geo)
    {
        fd=fopen(kmlfile,"wb");

        fprintf(fd,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        fprintf(fd,"<kml xmlns=\"http://earth.google.com/kml/2.1\">\n");
        fprintf(fd,"<!-- Generated by %s Version %s -->\n",SPLAT_NAME,SPLAT_VERSION);
        fprintf(fd,"  <Folder>\n");
        fprintf(fd,"   <name>%s</name>\n",SPLAT_NAME);
        fprintf(fd,"	 <description>%s Transmitter Contours</description>\n",xmtr[0].name);
        fprintf(fd,"	   <GroundOverlay>\n");
        fprintf(fd,"		 <name>SPLAT! Signal Power Level Contours</name>\n");
        fprintf(fd,"		 <description>SPLAT! Coverage</description>\n");
        fprintf(fd,"		 <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		 <Icon>\n");
        fprintf(fd,"			  <href>%s</href>\n",basename_s(mapfile));
        fprintf(fd,"		 </Icon>\n");
        /* fprintf(fd,"			<opacity>128</opacity>\n"); */
        fprintf(fd,"			<LatLonBox>\n");
        fprintf(fd,"			   <north>%.5f</north>\n",north);
        fprintf(fd,"			   <south>%.5f</south>\n",south);
        fprintf(fd,"			   <east>%.5f</east>\n",east);
        fprintf(fd,"			   <west>%.5f</west>\n",west);
        fprintf(fd,"			   <rotation>0.0</rotation>\n");
        fprintf(fd,"			</LatLonBox>\n");
        fprintf(fd,"	   </GroundOverlay>\n");
        fprintf(fd,"	   <ScreenOverlay>\n");
        fprintf(fd,"		  <name>Color Key</name>\n");
        fprintf(fd,"		  <description>Contour Color Key</description>\n");
        fprintf(fd,"		  <color>8cffffff</color>\n"); /* transparency level, higher is clearer */
        fprintf(fd,"		  <Icon>\n");
        fprintf(fd,"			<href>%s</href>\n",basename_s(ckfile));
        fprintf(fd,"		  </Icon>\n");
        fprintf(fd,"		  <overlayXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <screenXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <rotationXY x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"		  <size x=\"0\" y=\"0\" xunits=\"fraction\" yunits=\"fraction\"/>\n");
        fprintf(fd,"	   </ScreenOverlay>\n");

        for (unsigned int i=0; i<txsites; i++)
        {
            fprintf(fd,"	 <Placemark>\n");
            fprintf(fd,"	   <name>%s</name>\n",xmtr[i].name);
            fprintf(fd,"	   <visibility>1</visibility>\n");
            fprintf(fd,"	   <Style>\n");
            fprintf(fd,"	   <IconStyle>\n");
            fprintf(fd,"		<Icon>\n");
            fprintf(fd,"		  <href>root://icons/palette-5.png</href>\n");
            fprintf(fd,"		  <x>224</x>\n");
            fprintf(fd,"		  <y>224</y>\n");
            fprintf(fd,"		  <w>32</w>\n");
            fprintf(fd,"		  <h>32</h>\n");
            fprintf(fd,"		</Icon>\n");
            fprintf(fd,"	   </IconStyle>\n");
            fprintf(fd,"	   </Style>\n");
            fprintf(fd,"	  <Point>\n");
            fprintf(fd,"		<extrude>1</extrude>\n");
            fprintf(fd,"		<altitudeMode>relativeToGround</altitudeMode>\n");
            fprintf(fd,"		<coordinates>%f,%f,%f</coordinates>\n",(xmtr[i].lon<180.0?-xmtr[i].lon:360.0-xmtr[i].lon), xmtr[i].lat, xmtr[i].alt);
            fprintf(fd,"	  </Point>\n");
            fprintf(fd,"	 </Placemark>\n");
        }

        fprintf(fd,"  </Folder>\n");
        fprintf(fd,"</kml>\n");

        fclose(fd);
    }

    if (kml || geo)
    {
        /* No bottom legend */
        imgwidth = width;
        imgheight = height;
    }
    else
    {
        /* Allow space for bottom legend */
        imgwidth = width;
        imgheight = height + 30;
    }

    fprintf(stdout,"\nWriting Power Level (dBm) map \"%s\" (%ux%u image)... ",mapfile,imgwidth,imgheight);
    fflush(stdout);

    if (ImageWriterInit(&iw,mapfile,imagetype,imgwidth,imgheight)<0) {
        fprintf(stdout,"\nError writing \"%s\"!\n",mapfile);
        fflush(stdout);
        return;
    }

    for (y=0, lat=north; y<(int)height; y++, lat=north-(dpp*(double)y))
    {
        for (x=0, lon=max_west; x<(int)width; x++, lon=max_west-(dpp*(double)x))
        {
            if (lon<0.0)
                lon+=360.0;

            dem = FindDEM(lat, lon, x0, y0);
            if (dem)
            {
                mask=dem->mask[x0][y0];
                dBm=(dem->signal[x0][y0])-200;

                match=255;

                red=0;
                green=0;
                blue=0;

                if (dBm>=region.level[0])
                    match=0;
                else
                {
                    for (z=1; (z<region.levels && match==255); z++)
                    {
                        if (dBm<region.level[z-1] && dBm>=region.level[z])
                            match=z;
                    }
                }

                if (match<region.levels)
                {
                    if (smooth_contours && match>0)
                    {
                        red=(unsigned)interpolate(region.color[match][0],region.color[match-1][0],region.level[match],region.level[match-1],dBm);
                        green=(unsigned)interpolate(region.color[match][1],region.color[match-1][1],region.level[match],region.level[match-1],dBm);
                        blue=(unsigned)interpolate(region.color[match][2],region.color[match-1][2],region.level[match],region.level[match-1],dBm);
                    }

                    else
                    {
                        red=region.color[match][0];
                        green=region.color[match][1];
                        blue=region.color[match][2];
                    }
                }

                if (mask&2) 
                {
                    /* Text Labels: Red or otherwise */

                    if (red>=180 && green<=75 && blue<=75 && dBm!=0)
                        pixel=RGB(255^red,255^green,255^blue);
                    else
                        pixel=COLOR_RED;
                }
                else if (mask&4)
                {
                    /* County Boundaries: Black */
                    pixel=COLOR_BLACK;
                }
                else
                {
                    if (contour_threshold!=0 && dBm<contour_threshold)
                    {
                        if (ngs) /* No terrain */
                            pixel=COLOR_WHITE;
                        else
                        {
                            /* Display land or sea elevation */

                            if (dem->data[x0][y0]==0)
                                pixel=COLOR_MEDIUMBLUE;
                            else
                            {
                                terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                pixel=RGB(terrain,terrain,terrain);
                            }
                        }
                    }

                    else
                    {
                        /* Plot signal power level regions in color */

                        if (red!=0 || green!=0 || blue!=0)
                            pixel=RGB(red,green,blue);

                        else  /* terrain / sea-level */
                        {
                            if (ngs)
                                pixel=COLOR_WHITE;
                            else
                            {
                                if (dem->data[x0][y0]==0)
                                    pixel=COLOR_MEDIUMBLUE;
                                else
                                {
                                    /* Elevation: Greyscale */
                                    terrain=(unsigned)(0.5+pow((double)(dem->data[x0][y0]-min_elevation),one_over_gamma)*conversion);
                                    pixel=RGB(terrain,terrain,terrain);
                                }
                            }
                        }
                    }
                }
            }

            else
            {
                /* We should never get here, but if */
                /* we do, display the region as black */
                pixel=COLOR_BLACK;
            }

            ImageWriterAppendPixel(&iw, pixel);
        }

        ImageWriterEmitLine(&iw);
    }

    if (!kml && !geo)
    {
        /* Display legend along bottom of image
           if not generating .kml or .geo output. */

        colorwidth=(int)rint((float)width/(float)region.levels);

        for (y0=0; y0<30; y0++)
        {
            for (x0=0; x0<(int)width; x0++)
            {
                indx=x0/colorwidth;
                x=x0%colorwidth;

                level=abs(region.level[indx]);

                hundreds=level/100;

                if (hundreds>0)
                    level-=(hundreds*100);

                tens=level/10;

                if (tens>0)
                    level-=(tens*10);

                units=level;

                if (y0>=8 && y0<=23)
                {
                    if (hundreds>0)
                    {
                        if (region.level[indx]<0)
                        {
                            if (x>=5 && x<=12)
                                if (fontdata[16*('-')+(y0-8)]&(128>>(x-5)))
                                    indx=255;
                        }

                        else
                        {
                            if (x>=5 && x<=12)
                                if (fontdata[16*('+')+(y0-8)]&(128>>(x-5)))
                                    indx=255;
                        }

                        if (x>=13 && x<=20)	 
                            if (fontdata[16*(hundreds+'0')+(y0-8)]&(128>>(x-13)))
                                indx=255; 
                    }

                    if (tens>0 || hundreds>0)
                    {
                        if (hundreds==0)
                        {
                            if (region.level[indx]<0)
                            {
                                if (x>=13 && x<=20)
                                    if (fontdata[16*('-')+(y0-8)]&(128>>(x-13)))
                                        indx=255;
                            }

                            else
                            {
                                if (x>=13 && x<=20)
                                    if (fontdata[16*('+')+(y0-8)]&(128>>(x-13)))
                                        indx=255;
                            }
                        }

                        if (x>=21 && x<=28)	 
                            if (fontdata[16*(tens+'0')+(y0-8)]&(128>>(x-21)))
                                indx=255;
                    }

                    if (hundreds==0 && tens==0)
                    {
                        if (region.level[indx]<0)
                        {
                            if (x>=21 && x<=28)
                                if (fontdata[16*('-')+(y0-8)]&(128>>(x-21)))
                                    indx=255;
                        }

                        else
                        {
                            if (x>=21 && x<=28)
                                if (fontdata[16*('+')+(y0-8)]&(128>>(x-21)))
                                    indx=255;
                        }
                    }

                    if (x>=29 && x<=36)
                        if (fontdata[16*(units+'0')+(y0-8)]&(128>>(x-29)))
                            indx=255;

                    if (x>=37 && x<=44)
                        if (fontdata[16*('d')+(y0-8)]&(128>>(x-37)))
                            indx=255;

                    if (x>=45 && x<=52)
                        if (fontdata[16*('B')+(y0-8)]&(128>>(x-45)))
                            indx=255;

                    if (x>=53 && x<=60)
                        if (fontdata[16*('m')+(y0-8)]&(128>>(x-53)))
                            indx=255;
                }

                if (indx>region.levels)
                    pixel=COLOR_BLACK;
                else
                {
                    red=region.color[indx][0];
                    green=region.color[indx][1];
                    blue=region.color[indx][2];

                    pixel=RGB(red,green,blue);
                }

                ImageWriterAppendPixel(&iw, pixel);
            } 

            ImageWriterEmitLine(&iw);
        }
    }

    ImageWriterFinish(&iw);

    if (kml)
    {
        /* Write colorkey image file */
        height=30*region.levels;
        width=100;

        if (ImageWriterInit(&iw,ckfile,imagetype,width,height)<0) {
            fprintf(stdout,"\nError writing \"%s\"!\n",ckfile);
            fflush(stdout);
            return;
        }

        for (y0=0; y0<(int)height; y0++)
        {
            for (x0=0; x0<(int)width; x0++)
            {
                indx=y0/30;
                x=x0;

                level=abs(region.level[indx]);

                hundreds=level/100;

                if (hundreds>0)
                    level-=(hundreds*100);

                tens=level/10;

                if (tens>0)
                    level-=(tens*10);

                units=level;


                if ((y0%30)>=8 && (y0%30)<=23)
                {
                    if (hundreds>0)
                    {
                        if (region.level[indx]<0)
                        {
                            if (x>=5 && x<=12)
                                if (fontdata[16*('-')+((y0%30)-8)]&(128>>(x-5)))
                                    indx=255;
                        }

                        else
                        {
                            if (x>=5 && x<=12)
                                if (fontdata[16*('+')+((y0%30)-8)]&(128>>(x-5)))
                                    indx=255;
                        }

                        if (x>=13 && x<=20)	 
                            if (fontdata[16*(hundreds+'0')+((y0%30)-8)]&(128>>(x-13)))
                                indx=255; 
                    }

                    if (tens>0 || hundreds>0)
                    {
                        if (hundreds==0)
                        {
                            if (region.level[indx]<0)
                            {
                                if (x>=13 && x<=20)
                                    if (fontdata[16*('-')+((y0%30)-8)]&(128>>(x-13)))
                                        indx=255;
                            }

                            else
                            {
                                if (x>=13 && x<=20)
                                    if (fontdata[16*('+')+((y0%30)-8)]&(128>>(x-13)))
                                        indx=255;
                            }
                        }

                        if (x>=21 && x<=28)	 
                            if (fontdata[16*(tens+'0')+((y0%30)-8)]&(128>>(x-21)))
                                indx=255;
                    }

                    if (hundreds==0 && tens==0)
                    {
                        if (region.level[indx]<0)
                        {
                            if (x>=21 && x<=28)
                                if (fontdata[16*('-')+((y0%30)-8)]&(128>>(x-21)))
                                    indx=255;
                        }

                        else
                        {
                            if (x>=21 && x<=28)
                                if (fontdata[16*('+')+((y0%30)-8)]&(128>>(x-21)))
                                    indx=255;
                        }
                    }

                    if (x>=29 && x<=36)
                        if (fontdata[16*(units+'0')+((y0%30)-8)]&(128>>(x-29)))
                            indx=255;

                    if (x>=37 && x<=44)
                        if (fontdata[16*('d')+((y0%30)-8)]&(128>>(x-37)))
                            indx=255;

                    if (x>=45 && x<=52)
                        if (fontdata[16*('B')+((y0%30)-8)]&(128>>(x-45)))
                            indx=255;

                    if (x>=53 && x<=60)
                        if (fontdata[16*('m')+((y0%30)-8)]&(128>>(x-53)))
                            indx=255;
                }

                if (indx>region.levels) {
                    pixel=COLOR_BLACK;
                }
                else
                {
                    red=region.color[indx][0];
                    green=region.color[indx][1];
                    blue=region.color[indx][2];

                    pixel=RGB(red,green,blue);
                }

                ImageWriterAppendPixel(&iw, pixel);
            }

            ImageWriterEmitLine(&iw);
        }

        ImageWriterFinish(&iw);
    }

    free(mapfile);
    free(geofile);
    free(kmlfile);
    free(ckfile);

    fprintf(stdout,"Done!\n");
    fflush(stdout);
}

/* Invokes gnuplot to generate a file indicating the terrain profile between
 * the source and destination locations.
 * "basename" is the name assigned to the output file generated by gnuplot.
 * The filename extension is used to set gnuplot's terminal setting
 * and output file type.
 * If no extension is found, .png is assumed.
 */
void GraphTerrain(Site source, Site destination, char *name)
{
    int	x;
    char	basename[MAX_PATH_LEN], term[30], ext[20];
    double	minheight=100000.0, maxheight=-100000.0;
    FILE	*fd=NULL, *fd1=NULL;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return;
    }
    ReadPath(destination,source,path); 

    fd=fopen("profile.gp","wb");

    if (clutter>0.0)
        fd1=fopen("clutter.gp","wb");

    for (x=0; x<path->length; x++)
    {
        if ((path->elevation[x]+clutter)>maxheight)
            maxheight=path->elevation[x]+clutter;

        if (path->elevation[x]<minheight)
            minheight=path->elevation[x];

        if (metric)
        {
            fprintf(fd,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*path->elevation[x]);

            if (fd1!=NULL && x>0 && x<path->length-2)
                fprintf(fd1,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*(path->elevation[x]==0.0?path->elevation[x]:(path->elevation[x]+clutter)));
        }

        else
        {
            fprintf(fd,"%f\t%f\n",path->distance[x],path->elevation[x]);

            if (fd1!=NULL && x>0 && x<path->length-2)
                fprintf(fd1,"%f\t%f\n",path->distance[x],(path->elevation[x]==0.0?path->elevation[x]:(path->elevation[x]+clutter)));
        }
    }

    fclose(fd);

    if (fd1!=NULL)
        fclose(fd1);

    if (name[0]=='.')
    {
        /* Default filename and output file type */
        strcpy(basename,"profile");
        strcpy(ext,"png");
    }
    else
    {
        /* Extract extension and terminal type from "name" */
        strncpy(basename, name, MAX_PATH_LEN-1);
        stripExtension(basename, ext, 20);
        convertBackslashes(basename);
        if (strlen(ext) == 0)
        {
            strcpy(ext,"png");
        }
    }
    strcpy(term,ext);

    /* Either .ps or .postscript may be used
       as an extension for postscript output. */

    if (strncmp(term,"postscript",10)==0)
        strcpy(ext,"ps");

    else if (strncmp(ext,"ps",2)==0)
        strcpy(term,"postscript enhanced color");

    if (maxheight<1.0)
    {
        maxheight=1.0;	/* Avoid a gnuplot y-range error */ 
        minheight=-1.0;	/* over a completely sea-level path */
    }

    else
        minheight-=(0.01*maxheight);

    fd=fopen("splat.gp","w");
    fprintf(fd,"set grid\n");
    fprintf(fd,"set yrange [%2.3f to %2.3f]\n", metric?minheight*METERS_PER_FOOT:minheight, metric?maxheight*METERS_PER_FOOT:maxheight);
    fprintf(fd,"set encoding iso_8859_1\n");
    fprintf(fd,"set term %s\n",term);
    fprintf(fd,"set title \"%s Terrain Profile Between %s and %s (%.2f%c Azimuth)\"\n",SPLAT_NAME,destination.name, source.name, Azimuth(destination,source),176);

    if (metric)
    {
        fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f kilometers)\"\n",destination.name,source.name,KM_PER_MILE*Distance(source,destination));
        fprintf(fd,"set ylabel \"Ground Elevation Above Sea Level (meters)\"\n");
    }

    else
    {
        fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f miles)\"\n",destination.name,source.name,Distance(source,destination));
        fprintf(fd,"set ylabel \"Ground Elevation Above Sea Level (feet)\"\n");
    }

    fprintf(fd,"set output \"%s.%s\"\n",basename,ext);

    if (clutter>0.0)
    {
        if (metric)
            fprintf(fd,"plot \"profile.gp\" title \"Terrain Profile\" with lines, \"clutter.gp\" title \"Clutter Profile (%.2f meters)\" with lines\n",clutter*METERS_PER_FOOT);
        else
            fprintf(fd,"plot \"profile.gp\" title \"Terrain Profile\" with lines, \"clutter.gp\" title \"Clutter Profile (%.2f feet)\" with lines\n",clutter);
    }

    else
        fprintf(fd,"plot \"profile.gp\" title \"\" with lines\n");

    fclose(fd);

    x=system("gnuplot splat.gp");

    if (x!=-1)
    {
        if (!gpsav)
        {
            unlink("splat.gp");
            unlink("profile.gp");
        }

        fprintf(stdout,"Terrain plot written to: \"%s.%s\"\n",basename,ext);
        fflush(stdout);
    }

    else
        fprintf(stderr,"\n*** ERROR: Error occurred invoking gnuplot!\n");

    DestroyPath(path);
}

/* Invokes gnuplot to generate a file indicating the terrain elevation profile
 * between the source and destination locations.
 * "basename" is the name assigned to the output file generated by gnuplot.
 * The filename extension is used to set gnuplot's terminal setting
 * and output file type.
 * If no extension is found, .png is assumed.
 */
void GraphElevation(Site source, Site destination, char *name)
{
    int	x;
    char	basename[MAX_PATH_LEN], term[30], ext[20];
    double	angle, clutter_angle=0.0, refangle, maxangle=-90.0,
            minangle=90.0, distance;
    Site	remote, remote2;
    FILE	*fd=NULL, *fd1=NULL, *fd2=NULL;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return;
    }
    ReadPath(destination,source,path);  /* destination=RX, source=TX */
    refangle=ElevationAngle(destination,source);
    distance=Distance(source,destination);

    fd=fopen("profile.gp","wb");

    if (clutter>0.0)
        fd1=fopen("clutter.gp","wb");

    fd2=fopen("reference.gp","wb");

    for (x=1; x<path->length-1; x++)
    {
        remote.lat=path->lat[x];
        remote.lon=path->lon[x];
        remote.alt=0.0;
        angle=ElevationAngle(destination,remote);

        if (clutter>0.0)
        {
            remote2.lat=path->lat[x];
            remote2.lon=path->lon[x];

            if (path->elevation[x]!=0.0)
                remote2.alt=clutter;
            else
                remote2.alt=0.0;

            clutter_angle=ElevationAngle(destination,remote2);
        }

        if (metric)
        {
            fprintf(fd,"%f\t%f\n",KM_PER_MILE*path->distance[x],angle);

            if (fd1!=NULL)
                fprintf(fd1,"%f\t%f\n",KM_PER_MILE*path->distance[x],clutter_angle);

            fprintf(fd2,"%f\t%f\n",KM_PER_MILE*path->distance[x],refangle);
        }

        else
        {
            fprintf(fd,"%f\t%f\n",path->distance[x],angle);

            if (fd1!=NULL)
                fprintf(fd1,"%f\t%f\n",path->distance[x],clutter_angle);

            fprintf(fd2,"%f\t%f\n",path->distance[x],refangle);
        }

        if (angle>maxangle)
            maxangle=angle;

        if (clutter_angle>maxangle)
            maxangle=clutter_angle;

        if (angle<minangle)
            minangle=angle;
    }

    if (metric)
    {
        fprintf(fd,"%f\t%f\n",KM_PER_MILE*path->distance[path->length-1],refangle);
        fprintf(fd2,"%f\t%f\n",KM_PER_MILE*path->distance[path->length-1],refangle);
    }

    else
    {
        fprintf(fd,"%f\t%f\n",path->distance[path->length-1],refangle);
        fprintf(fd2,"%f\t%f\n",path->distance[path->length-1],refangle);
    }

    fclose(fd);

    if (fd1!=NULL)
        fclose(fd1);

    fclose(fd2);

    if (name[0]=='.')
    {
        /* Default filename and output file type */
        strcpy(basename,"profile");
        strcpy(ext,"png");
    }
    else
    {
        /* Extract extension and terminal type from "name" */
        strncpy(basename, name, MAX_PATH_LEN-1);
        stripExtension(basename, ext, 20);
        convertBackslashes(basename);
        if (strlen(ext) == 0)
        {
            strcpy(ext,"png");
        }
    }

    strcpy(term,ext);

    /* Either .ps or .postscript may be used
       as an extension for postscript output. */

    if (strncmp(term,"postscript",10)==0)
        strcpy(ext,"ps");

    else if (strncmp(ext,"ps",2)==0)
        strcpy(term,"postscript enhanced color");

    fd=fopen("splat.gp","w");

    fprintf(fd,"set grid\n");

    if (distance>2.0)
        fprintf(fd,"set yrange [%2.3f to %2.3f]\n", (-fabs(refangle)-0.25), maxangle+0.25);
    else
        fprintf(fd,"set yrange [%2.3f to %2.3f]\n", minangle, refangle+(-minangle/8.0));

    fprintf(fd,"set encoding iso_8859_1\n");
    fprintf(fd,"set term %s\n",term);
    fprintf(fd,"set title \"%s Elevation Profile Between %s and %s (%.2f%c azimuth)\"\n",SPLAT_NAME,destination.name,source.name,Azimuth(destination,source),176);

    if (metric)
        fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f kilometers)\"\n",destination.name,source.name,KM_PER_MILE*distance);
    else
        fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f miles)\"\n",destination.name,source.name,distance);


    fprintf(fd,"set ylabel \"Elevation Angle Along LOS Path Between\\n%s and %s (degrees)\"\n",destination.name,source.name);
    fprintf(fd,"set output \"%s.%s\"\n",basename,ext);

    if (clutter>0.0)
    {
        if (metric)
            fprintf(fd,"plot \"profile.gp\" title \"Real Earth Profile\" with lines, \"clutter.gp\" title \"Clutter Profile (%.2f meters)\" with lines, \"reference.gp\" title \"Line of Sight Path (%.2f%c elevation)\" with lines\n",clutter*METERS_PER_FOOT,refangle,176);
        else
            fprintf(fd,"plot \"profile.gp\" title \"Real Earth Profile\" with lines, \"clutter.gp\" title \"Clutter Profile (%.2f feet)\" with lines, \"reference.gp\" title \"Line of Sight Path (%.2f%c elevation)\" with lines\n",clutter,refangle,176);
    }

    else
        fprintf(fd,"plot \"profile.gp\" title \"Real Earth Profile\" with lines, \"reference.gp\" title \"Line of Sight Path (%.2f%c elevation)\" with lines\n",refangle,176);

    fclose(fd);

    x=system("gnuplot splat.gp");

    if (x!=-1)
    {
        if (!gpsav)
        {
            unlink("splat.gp");
            unlink("profile.gp");
            unlink("reference.gp");

            if (clutter>0.0)
                unlink("clutter.gp");
        }

        fprintf(stdout,"Elevation plot written to: \"%s.%s\"\n",basename,ext);
        fflush(stdout);
    }

    else
        fprintf(stderr,"\n*** ERROR: Error occurred invoking gnuplot!\n");

    DestroyPath(path);
}

/* Invokes gnuplot to generate a file indicating the terrain height profile
 * between the source and destination locations referenced to the lint-of-sight
 * path between the receive and transmit sites.
 * "basename" is the name assigned to the output file generated by gnuplot.
 * The filename extension is used to set gnuplot's terminal setting
 * and output file type.
 * If no extension is found, .png is assumed.
 */
void GraphHeight(Site source, Site destination, char *name, bool fresnel_plot, bool normalized)
{
    int	x;
    char	basename[MAX_PATH_LEN], term[30], ext[20];
    double	a, b, c, height=0.0, refangle, cangle, maxheight=-100000.0,
            minheight=100000.0, lambda=0.0, f_zone=0.0, fpt6_zone=0.0,
            nm=0.0, nb=0.0, ed=0.0, es=0.0, r=0.0, d=0.0, d1=0.0,
            terrain, azimuth, distance, dheight=0.0, minterrain=100000.0,
            minearth=100000.0, miny, maxy, min2y, max2y;
    Site	remote;
    FILE	*fd=NULL, *fd1=NULL, *fd2=NULL, *fd3=NULL, *fd4=NULL, *fd5=NULL;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return;
    }
    ReadPath(destination,source,path);  /* destination=RX, source=TX */
    azimuth=Azimuth(destination,source);
    distance=Distance(destination,source);
    refangle=ElevationAngle(destination,source);

    if (destination.amsl_flag)
        b=destination.alt+earthradius;
    else
        b=GetElevation(destination)+destination.alt+earthradius;

    /* Wavelength and path distance (great circle) in feet. */

    if (fresnel_plot)
    {
        lambda=9.8425e8/(LR.frq_mhz*1e6);
        d=5280.0*path->distance[path->length-1];
    }

    if (normalized)
    {
        ed=GetElevation(destination);
        es=GetElevation(source);
        nb=-destination.alt-ed;
        nm=(-source.alt-es-nb)/(path->distance[path->length-1]);
    }

    fd=fopen("profile.gp","wb");

    if (clutter>0.0)
        fd1=fopen("clutter.gp","wb");

    fd2=fopen("reference.gp","wb");
    fd5=fopen("curvature.gp", "wb");

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
    {
        fd3=fopen("fresnel.gp", "wb");
        fd4=fopen("fresnel_pt_6.gp", "wb");
    }

    for (x=0; x<path->length-1; x++)
    {
        remote.lat=path->lat[x];
        remote.lon=path->lon[x];
        remote.alt=0.0;

        terrain=GetElevation(remote);

        if (x==0)
            terrain+=destination.alt;  /* RX antenna spike */

        a=terrain+earthradius;
        cangle=5280.0*Distance(destination,remote)/earthradius;
        c=b*sin(refangle*DEG2RAD+HALFPI)/sin(HALFPI-refangle*DEG2RAD-cangle);

        height=a-c;

        /* Per Fink and Christiansen, Electronics
         * Engineers' Handbook, 1989:
         *
         *   H = sqrt(lamba * d1 * (d - d1)/d)
         *
         * where H is the distance from the LOS
         * path to the first Fresnel zone boundary.
         */

        if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
        {
            d1=5280.0*path->distance[x];
            f_zone=-1.0*sqrt(lambda*d1*(d-d1)/d);
            fpt6_zone=f_zone*fzone_clearance;
        }

        if (normalized)
        {
            r=-(nm*path->distance[x])-nb;
            height+=r;

            if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
            {
                f_zone+=r;
                fpt6_zone+=r;
            }
        }

        else
            r=0.0;

        if (metric)
        {
            fprintf(fd,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*height);

            if (fd1!=NULL && x>0 && x<path->length-2)
                fprintf(fd1,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*(terrain==0.0?height:(height+clutter)));

            fprintf(fd2,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*r);
            fprintf(fd5,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*(height-terrain));
        }

        else
        {
            fprintf(fd,"%f\t%f\n",path->distance[x],height);

            if (fd1!=NULL && x>0 && x<path->length-2)
                fprintf(fd1,"%f\t%f\n",path->distance[x],(terrain==0.0?height:(height+clutter)));

            fprintf(fd2,"%f\t%f\n",path->distance[x],r);
            fprintf(fd5,"%f\t%f\n",path->distance[x],height-terrain);
        }

        if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
        {
            if (metric)
            {
                fprintf(fd3,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*f_zone);
                fprintf(fd4,"%f\t%f\n",KM_PER_MILE*path->distance[x],METERS_PER_FOOT*fpt6_zone);
            }

            else
            {
                fprintf(fd3,"%f\t%f\n",path->distance[x],f_zone);
                fprintf(fd4,"%f\t%f\n",path->distance[x],fpt6_zone);
            }

            if (f_zone<minheight)
                minheight=f_zone;
        }

        if ((height+clutter)>maxheight)
            maxheight=height+clutter;

        if (height<minheight)
            minheight=height;

        if (r>maxheight)
            maxheight=r;

        if (terrain<minterrain)
            minterrain=terrain;

        if ((height-terrain)<minearth)
            minearth=height-terrain;
    }

    if (normalized)
        r=-(nm*path->distance[path->length-1])-nb;
    else
        r=0.0;

    if (metric)
    {
        fprintf(fd,"%f\t%f\n",KM_PER_MILE*path->distance[path->length-1],METERS_PER_FOOT*r);
        fprintf(fd2,"%f\t%f\n",KM_PER_MILE*path->distance[path->length-1],METERS_PER_FOOT*r);
    }

    else
    {
        fprintf(fd,"%f\t%f\n",path->distance[path->length-1],r);
        fprintf(fd2,"%f\t%f\n",path->distance[path->length-1],r);
    }

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
    {
        if (metric)
        {
            fprintf(fd3,"%f\t%f\n",KM_PER_MILE*path->distance[path->length-1],METERS_PER_FOOT*r);
            fprintf(fd4,"%f\t%f\n",KM_PER_MILE*path->distance[path->length-1],METERS_PER_FOOT*r);
        }

        else
        {
            fprintf(fd3,"%f\t%f\n",path->distance[path->length-1],r);
            fprintf(fd4,"%f\t%f\n",path->distance[path->length-1],r);
        }
    }

    if (r>maxheight)
        maxheight=r;

    if (r<minheight)
        minheight=r;

    fclose(fd);

    if (fd1!=NULL)
        fclose(fd1);

    fclose(fd2);
    fclose(fd5);

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
    {
        fclose(fd3);
        fclose(fd4);
    }

    if (name[0]=='.')
    {
        /* Default filename and output file type */
        strcpy(basename,"profile");
        strcpy(ext,"png");
    }
    else
    {
        /* Extract extension and terminal type from "name" */
        strncpy(basename, name, MAX_PATH_LEN-1);
        stripExtension(basename, ext, 20);
        convertBackslashes(basename);
        if (strlen(ext) == 0)
        {
            strcpy(ext,"png");
        }
    }
    strcpy(term,ext);

    /* Either .ps or .postscript may be used
       as an extension for postscript output. */

    if (strncmp(term,"postscript",10)==0)
        strcpy(ext,"ps");

    else if (strncmp(ext,"ps",2)==0)
        strcpy(term,"postscript enhanced color");

    fd=fopen("splat.gp","w");

    dheight=maxheight-minheight;
    miny=minheight-0.15*dheight;
    maxy=maxheight+0.05*dheight;

    if (maxy<20.0)
        maxy=20.0;

    dheight=maxheight-minheight;
    min2y=miny-minterrain+0.05*dheight;

    if (minearth<min2y)
    {
        miny-=min2y-minearth+0.05*dheight;
        min2y=minearth-0.05*dheight;
    }

    max2y=min2y+maxy-miny;

    fprintf(fd,"set grid\n");
    fprintf(fd,"set yrange [%2.3f to %2.3f]\n", metric?miny*METERS_PER_FOOT:miny, metric?maxy*METERS_PER_FOOT:maxy);
    fprintf(fd,"set y2range [%2.3f to %2.3f]\n", metric?min2y*METERS_PER_FOOT:min2y, metric?max2y*METERS_PER_FOOT:max2y);
    fprintf(fd,"set xrange [-0.5 to %2.3f]\n",metric?KM_PER_MILE*rint(distance+0.5):rint(distance+0.5));
    fprintf(fd,"set encoding iso_8859_1\n");
    fprintf(fd,"set term %s\n",term);

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
        fprintf(fd,"set title \"%s Path Profile Between %s and %s (%.2f%c azimuth)\\nWith First Fresnel Zone\"\n",SPLAT_NAME, destination.name, source.name, azimuth,176);

    else
        fprintf(fd,"set title \"%s Height Profile Between %s and %s (%.2f%c azimuth)\"\n",SPLAT_NAME, destination.name, source.name, azimuth,176);

    if (metric)
        fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f kilometers)\"\n",destination.name,source.name,KM_PER_MILE*Distance(source,destination));
    else
        fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f miles)\"\n",destination.name,source.name,Distance(source,destination));

    if (normalized)
    {
        if (metric)
            fprintf(fd,"set ylabel \"Normalized Height Referenced To LOS Path Between\\n%s and %s (meters)\"\n",destination.name,source.name);

        else
            fprintf(fd,"set ylabel \"Normalized Height Referenced To LOS Path Between\\n%s and %s (feet)\"\n",destination.name,source.name);

    }

    else
    {
        if (metric)
            fprintf(fd,"set ylabel \"Height Referenced To LOS Path Between\\n%s and %s (meters)\"\n",destination.name,source.name);

        else
            fprintf(fd,"set ylabel \"Height Referenced To LOS Path Between\\n%s and %s (feet)\"\n",destination.name,source.name);
    }

    fprintf(fd,"set output \"%s.%s\"\n",basename,ext);

    if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
    {
        if (clutter>0.0)
        {
            if (metric)
                fprintf(fd,"plot \"profile.gp\" title \"Point-to-Point Profile\" with lines, \"clutter.gp\" title \"Ground Clutter (%.2f meters)\" with lines, \"reference.gp\" title \"Line of Sight Path\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature Contour\" with lines, \"fresnel.gp\" axes x1y1 title \"First Fresnel Zone (%.3f MHz)\" with lines, \"fresnel_pt_6.gp\" title \"%.0f%% of First Fresnel Zone\" with lines\n",clutter*METERS_PER_FOOT,LR.frq_mhz,fzone_clearance*100.0);
            else
                fprintf(fd,"plot \"profile.gp\" title \"Point-to-Point Profile\" with lines, \"clutter.gp\" title \"Ground Clutter (%.2f feet)\" with lines, \"reference.gp\" title \"Line of Sight Path\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature Contour\" with lines, \"fresnel.gp\" axes x1y1 title \"First Fresnel Zone (%.3f MHz)\" with lines, \"fresnel_pt_6.gp\" title \"%.0f%% of First Fresnel Zone\" with lines\n",clutter,LR.frq_mhz,fzone_clearance*100.0);
        }

        else
            fprintf(fd,"plot \"profile.gp\" title \"Point-to-Point Profile\" with lines, \"reference.gp\" title \"Line of Sight Path\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature Contour\" with lines, \"fresnel.gp\" axes x1y1 title \"First Fresnel Zone (%.3f MHz)\" with lines, \"fresnel_pt_6.gp\" title \"%.0f%% of First Fresnel Zone\" with lines\n",LR.frq_mhz,fzone_clearance*100.0);
    }

    else
    {
        if (clutter>0.0)
        {
            if (metric)
                fprintf(fd,"plot \"profile.gp\" title \"Point-to-Point Profile\" with lines, \"clutter.gp\" title \"Ground Clutter (%.2f meters)\" with lines, \"reference.gp\" title \"Line Of Sight Path\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature Contour\" with lines\n",clutter*METERS_PER_FOOT);
            else
                fprintf(fd,"plot \"profile.gp\" title \"Point-to-Point Profile\" with lines, \"clutter.gp\" title \"Ground Clutter (%.2f feet)\" with lines, \"reference.gp\" title \"Line Of Sight Path\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature Contour\" with lines\n",clutter);
        }

        else
            fprintf(fd,"plot \"profile.gp\" title \"Point-to-Point Profile\" with lines, \"reference.gp\" title \"Line Of Sight Path\" with lines, \"curvature.gp\" axes x1y2 title \"Earth's Curvature Contour\" with lines\n");

    }

    fclose(fd);

    x=system("gnuplot splat.gp");

    if (x!=-1)
    {
        if (!gpsav)
        {
            unlink("splat.gp");
            unlink("profile.gp");
            unlink("reference.gp");
            unlink("curvature.gp");

            if (fd1!=NULL)
                unlink("clutter.gp");

            if ((LR.frq_mhz>=20.0) && (LR.frq_mhz<=20000.0) && fresnel_plot)
            {
                unlink("fresnel.gp");
                unlink("fresnel_pt_6.gp");
            }
        }

        fprintf(stdout,"\n%seight plot written to: \"%s.%s\"",
                normalized?"Normalized h":"H",basename,ext);
        fflush(stdout);
    }

    else
        fprintf(stderr,"\n*** ERROR: Error occurred invoking gnuplot!\n");

    DestroyPath(path);
}

/* Performs an obstruction analysis along the path between receiver
 * and transmitter.
 */
void ObstructionAnalysis(Site xmtr, Site rcvr, double f, FILE *outfile)
{
    int	x;
    Site site_x;
    double	h_r, h_t, h_x, h_r_orig, cos_tx_angle, cos_test_angle,
            cos_tx_angle_f1, cos_tx_angle_fpt6, d_tx, d_x,
            h_r_f1, h_r_fpt6, h_f, h_los, lambda=0.0;
    char	buf[MAX_PATH_LEN], buf_fpt6[MAX_PATH_LEN], buf_f1[MAX_PATH_LEN];

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return;
    }
    ReadPath(xmtr,rcvr,path);

    if (rcvr.amsl_flag)
        h_r=rcvr.alt+earthradius;
    else
        h_r=GetElevation(rcvr)+rcvr.alt+earthradius;

    h_r_f1=h_r;
    h_r_fpt6=h_r;
    h_r_orig=h_r;

    if (xmtr.amsl_flag)
        h_t=xmtr.alt+earthradius;
    else
        h_t=GetElevation(xmtr)+xmtr.alt+earthradius;

    d_tx=5280.0*Distance(rcvr,xmtr);
    cos_tx_angle=((h_r*h_r)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r*d_tx);
    cos_tx_angle_f1=cos_tx_angle;
    cos_tx_angle_fpt6=cos_tx_angle;

    if (f)
        lambda=9.8425e8/(f*1e6);

    if (clutter>0.0)
    {
        fprintf(outfile,"Terrain has been raised by");

        if (metric)
            fprintf(outfile," %.2f meters",METERS_PER_FOOT*clutter);
        else
            fprintf(outfile," %.2f feet",clutter);

        fprintf(outfile," to account for ground clutter.\n\n");
    }

    /* At each point along the path calculate the cosine
       of a sort of "inverse elevation angle" at the receiver.
       From the antenna, 0 deg. looks at the ground, and 90 deg.
       is parallel to the ground.

       Start at the receiver.  If this is the lowest antenna,
       then terrain obstructions will be nearest to it.  (Plus,
       that's the way SPLAT!'s original los() did it.)

       Calculate cosines only.  That's sufficient to compare
       angles and it saves the extra computational burden of
       acos().  However, note the inverted comparison: if
       acos(A) > acos(B), then B > A. */

    for (x=path->length-1; x>0; x--)
    {
        site_x.lat=path->lat[x];
        site_x.lon=path->lon[x];
        site_x.alt=0.0;

        h_x=GetElevation(site_x)+earthradius+clutter;
        d_x=5280.0*Distance(rcvr,site_x);

        /* Deal with the LOS path first. */

        cos_test_angle=((h_r*h_r)+(d_x*d_x)-(h_x*h_x))/(2.0*h_r*d_x);

        if (cos_tx_angle>cos_test_angle)
        {
            if (h_r==h_r_orig)
                fprintf(outfile,"Between %s and %s, %s detected obstructions at:\n\n",rcvr.name,xmtr.name,SPLAT_NAME);

            if (site_x.lat>=0.0)
            {
                if (metric)
                    fprintf(outfile,"   %8.4f N,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",site_x.lat, site_x.lon, KM_PER_MILE*(d_x/5280.0), METERS_PER_FOOT*(h_x-earthradius));
                else
                    fprintf(outfile,"   %8.4f N,%9.4f W, %5.2f miles, %6.2f feet AMSL\n",site_x.lat, site_x.lon, d_x/5280.0, h_x-earthradius);
            }

            else
            {
                if (metric)
                    fprintf(outfile,"   %8.4f S,%9.4f W, %5.2f kilometers, %6.2f meters AMSL\n",-site_x.lat, site_x.lon, KM_PER_MILE*(d_x/5280.0), METERS_PER_FOOT*(h_x-earthradius));
                else

                    fprintf(outfile,"   %8.4f S,%9.4f W, %5.2f miles, %6.2f feet AMSL\n",-site_x.lat, site_x.lon, d_x/5280.0, h_x-earthradius);
            }
        }

        while (cos_tx_angle>cos_test_angle)
        {
            h_r+=1;
            cos_test_angle=((h_r*h_r)+(d_x*d_x)-(h_x*h_x))/(2.0*h_r*d_x);
            cos_tx_angle=((h_r*h_r)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r*d_tx);
        }

        if (f)
        {
            /* Now clear the first Fresnel zone... */

            cos_tx_angle_f1=((h_r_f1*h_r_f1)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_f1*d_tx);
            h_los=sqrt(h_r_f1*h_r_f1+d_x*d_x-2*h_r_f1*d_x*cos_tx_angle_f1);
            h_f=h_los-sqrt(lambda*d_x*(d_tx-d_x)/d_tx);

            while (h_f<h_x)
            {
                h_r_f1+=1;
                cos_tx_angle_f1=((h_r_f1*h_r_f1)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_f1*d_tx);
                h_los=sqrt(h_r_f1*h_r_f1+d_x*d_x-2*h_r_f1*d_x*cos_tx_angle_f1);
                h_f=h_los-sqrt(lambda*d_x*(d_tx-d_x)/d_tx);
            }

            /* and clear the 60% F1 zone. */

            cos_tx_angle_fpt6=((h_r_fpt6*h_r_fpt6)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_fpt6*d_tx);
            h_los=sqrt(h_r_fpt6*h_r_fpt6+d_x*d_x-2*h_r_fpt6*d_x*cos_tx_angle_fpt6);
            h_f=h_los-fzone_clearance*sqrt(lambda*d_x*(d_tx-d_x)/d_tx);

            while (h_f<h_x)
            {
                h_r_fpt6+=1;
                cos_tx_angle_fpt6=((h_r_fpt6*h_r_fpt6)+(d_tx*d_tx)-(h_t*h_t))/(2.0*h_r_fpt6*d_tx);
                h_los=sqrt(h_r_fpt6*h_r_fpt6+d_x*d_x-2*h_r_fpt6*d_x*cos_tx_angle_fpt6);
                h_f=h_los-fzone_clearance*sqrt(lambda*d_x*(d_tx-d_x)/d_tx);
            }
        }
    }

    if (h_r>h_r_orig)
    {
        if (metric)
            snprintf(buf,MAX_PATH_LEN,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear all obstructions detected by %s.\n",rcvr.name, METERS_PER_FOOT*(h_r-GetElevation(rcvr)-earthradius),SPLAT_NAME);
        else
            snprintf(buf,MAX_PATH_LEN,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear all obstructions detected by %s.\n",rcvr.name, h_r-GetElevation(rcvr)-earthradius,SPLAT_NAME);
    }

    else
        snprintf(buf,MAX_PATH_LEN,"\nNo obstructions to LOS path due to terrain were detected by %s\n",SPLAT_NAME);

    if (f)
    {
        if (h_r_fpt6>h_r_orig)
        {
            if (metric)
                snprintf(buf_fpt6,MAX_PATH_LEN,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear %.0f%c of the first Fresnel zone.\n",rcvr.name, METERS_PER_FOOT*(h_r_fpt6-GetElevation(rcvr)-earthradius),fzone_clearance*100.0,37);

            else
                snprintf(buf_fpt6,MAX_PATH_LEN,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear %.0f%c of the first Fresnel zone.\n",rcvr.name, h_r_fpt6-GetElevation(rcvr)-earthradius,fzone_clearance*100.0,37);
        }

        else
            snprintf(buf_fpt6,MAX_PATH_LEN,"\n%.0f%c of the first Fresnel zone is clear.\n",fzone_clearance*100.0,37);

        if (h_r_f1>h_r_orig)
        {
            if (metric)
                snprintf(buf_f1,MAX_PATH_LEN,"\nAntenna at %s must be raised to at least %.2f meters AGL\nto clear the first Fresnel zone.\n",rcvr.name, METERS_PER_FOOT*(h_r_f1-GetElevation(rcvr)-earthradius));

            else
                snprintf(buf_f1,MAX_PATH_LEN,"\nAntenna at %s must be raised to at least %.2f feet AGL\nto clear the first Fresnel zone.\n",rcvr.name, h_r_f1-GetElevation(rcvr)-earthradius);

        }

        else
            snprintf(buf_f1,MAX_PATH_LEN,"\nThe first Fresnel zone is clear.\n");
    }

    fprintf(outfile,"%s",buf);

    if (f)
    {
        fprintf(outfile,"%s",buf_f1);
        fprintf(outfile,"%s",buf_fpt6);
    }

    DestroyPath(path);
}

/* Writes a SPLAT! Path Report (source.name-to-dest.name.txt) to the filesystem.
 * If (graph_it == true), then gnuplot is invoked to generate an appropriate output
 * file indicating the ITM/ITWOM model loss between the source and destination
 * locations. You don't get a choice in the name of the text report, but "name"
 * is the name assigned to the output file generated by gnuplot.
 * The filename extension is used to set gnuplot's terminal setting and output
 * file type. If no extension is found, .png is assumed.
 */
void PathReport(Site source, Site destination, char *name, bool graph_it)
{
    int	x, y, errnum;
    char	basename[MAX_PATH_LEN], term[30], ext[20], strmode[100],
            report_name[MAX_PATH_LEN], block=0, propbuf[20];
    double	maxloss=-100000.0, minloss=100000.0, loss, haavt,
            angle1, angle2, azimuth, pattern=1.0, patterndB=0.0,
            total_loss=0.0, cos_xmtr_angle, cos_test_angle=0.0,
            source_alt, test_alt, dest_alt, source_alt2, dest_alt2,
            distance, elevation, four_thirds_earth, field_strength,
            free_space_loss=0.0, eirp=0.0, voltage, rxp, dBm,
            power_density;
    FILE	*fd=NULL, *fd2=NULL;

    snprintf(report_name,MAX_PATH_LEN,"%s-to-%s.txt",source.name,destination.name);

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path!\n");
        return;
    }

    four_thirds_earth=FOUR_THIRDS*EARTHRADIUS;

    for (x=0; report_name[x]!=0; x++)
        if (report_name[x]==32 || report_name[x]==17 || report_name[x]==92 || report_name[x]==42 || report_name[x]==47)
            report_name[x]='_';

    fd2=fopen(report_name,"w");

    fprintf(fd2,"\n\t\t--==[ %s v%s Path Analysis ]==--\n\n",SPLAT_NAME,SPLAT_VERSION);
    fprintf(fd2,"%s\n\n",dashes);
    fprintf(fd2,"Transmitter site: %s\n",source.name);

    if (source.lat>=0.0)
    {
        fprintf(fd2,"Site location: %.4f North / %.4f West",source.lat, source.lon);
        fprintf(fd2, " (%s N / ", dec2dms(source.lat));
    }

    else
    {

        fprintf(fd2,"Site location: %.4f South / %.4f West",-source.lat, source.lon);
        fprintf(fd2, " (%s S / ", dec2dms(source.lat));
    }

    fprintf(fd2, "%s W)\n", dec2dms(source.lon));

    if (metric)
    {
        if (source.amsl_flag)
            fprintf(fd2,"Antenna height: %.2f meters AMSL\n",METERS_PER_FOOT*source.alt);
        else
        {
            fprintf(fd2,"Ground elevation: %.2f meters AMSL\n",METERS_PER_FOOT*GetElevation(source));
            fprintf(fd2,"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",METERS_PER_FOOT*source.alt,METERS_PER_FOOT*(source.alt+GetElevation(source)));
        }
    }
    else
    {
        if (source.amsl_flag)
            fprintf(fd2,"Antenna height: %.2f feet AMSL\n",source.alt);
        else
        {
            fprintf(fd2,"Ground elevation: %.2f feet AMSL\n",GetElevation(source));
            fprintf(fd2,"Antenna height: %.2f feet AGL / %.2f feet AMSL\n",source.alt, source.alt+GetElevation(source));
        }
    }

    haavt=haat(source);

    if (haavt>-4999.0)
    {
        if (metric)
            fprintf(fd2,"Antenna height above average terrain: %.2f meters\n",METERS_PER_FOOT*haavt);
        else
            fprintf(fd2,"Antenna height above average terrain: %.2f feet\n",haavt);
    }

    azimuth=Azimuth(source,destination);
    angle1=ElevationAngle(source,destination);
    angle2=ElevationAngle2(source,destination,earthradius);

    if (got_azimuth_pattern || got_elevation_pattern)
    {
        x=(int)rint(10.0*(10.0-angle2));

        if (x>=0 && x<=1000) {
            pattern=(double)LR.antenna_pattern[(int)rint(azimuth)][x];

            if (pattern > 0.0) {
                patterndB=20.0*log10(pattern);
            } else if (pattern != NO_ANTENNA_DATA) {
                patterndB=-9999;
            }
        }
    }

    if (metric)
        fprintf(fd2,"Distance to %s: %.2f kilometers\n",destination.name,KM_PER_MILE*Distance(source,destination));

    else
        fprintf(fd2,"Distance to %s: %.2f miles\n",destination.name,Distance(source,destination));

    fprintf(fd2,"Azimuth to %s: %.2f degrees\n",destination.name,azimuth);

    if (angle1>=0.0)
        fprintf(fd2,"Elevation angle to %s: %+.4f degrees\n",destination.name,angle1);

    else
        fprintf(fd2,"Depression angle to %s: %+.4f degrees\n",destination.name,angle1);

    if ((angle2-angle1)>0.0001)
    {
        if (angle2<0.0)
            fprintf(fd2,"Depression");
        else
            fprintf(fd2,"Elevation");

        fprintf(fd2," angle to the first obstruction: %+.4f degrees\n",angle2);
    }

    fprintf(fd2,"\n%s\n\n",dashes);

    /* Receiver */

    fprintf(fd2,"Receiver site: %s\n",destination.name);

    if (destination.lat>=0.0)
    {
        fprintf(fd2,"Site location: %.4f North / %.4f West",destination.lat, destination.lon);
        fprintf(fd2, " (%s N / ", dec2dms(destination.lat));
    }

    else
    {
        fprintf(fd2,"Site location: %.4f South / %.4f West",-destination.lat, destination.lon);
        fprintf(fd2, " (%s S / ", dec2dms(destination.lat));
    }

    fprintf(fd2, "%s W)\n", dec2dms(destination.lon));

    if (metric)
    {
        if (destination.amsl_flag)
            fprintf(fd2,"Antenna height: %.2f meters AMSL\n", METERS_PER_FOOT*destination.alt);

        else
        {
            fprintf(fd2,"Ground elevation: %.2f meters AMSL\n",METERS_PER_FOOT*GetElevation(destination));
            fprintf(fd2,"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",METERS_PER_FOOT*destination.alt, METERS_PER_FOOT*(destination.alt+GetElevation(destination)));
        }
    }
    else
    {
        if (destination.amsl_flag)
            fprintf(fd2,"Antenna height: %.2f feet AMSL\n", destination.alt);
        else
        {
            fprintf(fd2,"Ground elevation: %.2f feet AMSL\n",GetElevation(destination));
            fprintf(fd2,"Antenna height: %.2f feet AGL / %.2f feet AMSL\n",destination.alt, destination.alt+GetElevation(destination));
        }
    }

    haavt=haat(destination);

    if (haavt>-4999.0)
    {
        if (metric)
            fprintf(fd2,"Antenna height above average terrain: %.2f meters\n",METERS_PER_FOOT*haavt);
        else
            fprintf(fd2,"Antenna height above average terrain: %.2f feet\n",haavt);
    }

    if (metric)
        fprintf(fd2,"Distance to %s: %.2f kilometers\n",source.name,KM_PER_MILE*Distance(source,destination));

    else
        fprintf(fd2,"Distance to %s: %.2f miles\n",source.name,Distance(source,destination));

    azimuth=Azimuth(destination,source);

    angle1=ElevationAngle(destination,source);
    angle2=ElevationAngle2(destination,source,earthradius);

    fprintf(fd2,"Azimuth to %s: %.2f degrees\n",source.name,azimuth);

    if (angle1>=0.0)
        fprintf(fd2,"Elevation angle to %s: %+.4f degrees\n",source.name,angle1);

    else
        fprintf(fd2,"Depression angle to %s: %+.4f degrees\n",source.name,angle1);

    if ((angle2-angle1)>0.0001)
    {
        if (angle2<0.0)
            fprintf(fd2,"Depression");
        else
            fprintf(fd2,"Elevation");

        fprintf(fd2," angle to the first obstruction: %+.4f degrees\n",angle2);
    }

    fprintf(fd2,"\n%s\n\n",dashes);

    if (LR.frq_mhz>0.0)
    {
        if (itwom)
            fprintf(fd2,"ITWOM Version %.1f Parameters Used In This Analysis:\n\n",ITWOMVersion());
        else
            fprintf(fd2,"Longley-Rice Parameters Used In This Analysis:\n\n");

        fprintf(fd2,"Earth's Dielectric Constant: %.3lf\n",LR.eps_dielect);
        fprintf(fd2,"Earth's Conductivity: %.3lf Siemens/meter\n",LR.sgm_conductivity);
        fprintf(fd2,"Atmospheric Bending Constant (N-units): %.3lf ppm\n",LR.eno_ns_surfref);
        fprintf(fd2,"Frequency: %.3lf MHz\n",LR.frq_mhz);
        fprintf(fd2,"Radio Climate: %d (",LR.radio_climate);

        switch (LR.radio_climate)
        {
            case 1:
                fprintf(fd2,"Equatorial");
                break;

            case 2:
                fprintf(fd2,"Continental Subtropical");
                break;

            case 3:
                fprintf(fd2,"Maritime Subtropical");
                break;

            case 4:
                fprintf(fd2,"Desert");
                break;

            case 5:
                fprintf(fd2,"Continental Temperate");
                break;

            case 6:
                fprintf(fd2,"Martitime Temperate, Over Land");
                break;

            case 7:
                fprintf(fd2,"Maritime Temperate, Over Sea");
                break;

            default:
                fprintf(fd2,"Unknown");
        }

        fprintf(fd2,")\nPolarization: %d (",LR.pol);

        if (LR.pol==0)
            fprintf(fd2,"Horizontal");

        if (LR.pol==1)
            fprintf(fd2,"Vertical");

        fprintf(fd2,")\nFraction of Situations: %.1lf%c\n",LR.conf*100.0,37);
        fprintf(fd2,"Fraction of Time: %.1lf%c\n",LR.rel*100.0,37);

        if (LR.erp!=0.0)
        {
            fprintf(fd2,"Transmitter ERP: ");

            if (LR.erp<1.0)
                fprintf(fd2,"%.1lf milliwatts",1000.0*LR.erp);

            if (LR.erp>=1.0 && LR.erp<10.0)
                fprintf(fd2,"%.1lf Watts",LR.erp);

            if (LR.erp>=10.0 && LR.erp<10.0e3)
                fprintf(fd2,"%.0lf Watts",LR.erp);

            if (LR.erp>=10.0e3)
                fprintf(fd2,"%.3lf kilowatts",LR.erp/1.0e3);

            dBm=10.0*(log10(LR.erp*1000.0));
            fprintf(fd2," (%+.2f dBm)\n",dBm);

            /* EIRP = ERP + 2.14 dB */

            fprintf(fd2,"Transmitter EIRP: ");

            eirp=LR.erp*1.636816521;

            if (eirp<1.0)
                fprintf(fd2,"%.1lf milliwatts",1000.0*eirp);

            if (eirp>=1.0 && eirp<10.0)
                fprintf(fd2,"%.1lf Watts",eirp);

            if (eirp>=10.0 && eirp<10.0e3)
                fprintf(fd2,"%.0lf Watts",eirp);

            if (eirp>=10.0e3)
                fprintf(fd2,"%.3lf kilowatts",eirp/1.0e3);

            dBm=10.0*(log10(eirp*1000.0));
            fprintf(fd2," (%+.2f dBm)\n",dBm);
        }

        fprintf(fd2,"\n%s\n\n",dashes);

        fprintf(fd2,"Summary For The Link Between %s and %s:\n\n",source.name, destination.name);

        if (patterndB!=0.0) {
            if (patterndB == -9999) {
                fprintf(fd2,"%s antenna pattern towards %s: -infinity (blocked)\n", source.name, destination.name);
            } else {
                fprintf(fd2,"%s antenna pattern towards %s: %.3f (%.2f dB)\n", source.name, destination.name, pattern, patterndB);
            }
        }

        ReadPath(source,destination,path);  /* source=TX, destination=RX */

        elev_t *elev = (elev_t*)calloc(path->arysize + 10, sizeof(elev_t));

        /* Copy elevations plus clutter along
           path into the elev[] array. */

        for (x=1; x<path->length-1; x++)
            elev[x+2]=METERS_PER_FOOT*(path->elevation[x]==0.0?(elev_t)(path->elevation[x]):(elev_t)(clutter+path->elevation[x]));

        /* Copy ending points without clutter */

        elev[2]=(elev_t)(path->elevation[0]*METERS_PER_FOOT);
        elev[path->length+1]=(elev_t)(path->elevation[path->length-1]*METERS_PER_FOOT);

        fd=fopen("profile.gp","w");

        azimuth=rint(Azimuth(source,destination));

        for (y=2; y<(path->length-1); y++)  /* path->length-1 avoids LR error */
        {
            distance=5280.0*path->distance[y];

            if (source.amsl_flag)
                source_alt=four_thirds_earth+source.alt;
            else
                source_alt=four_thirds_earth+source.alt+path->elevation[0];

            /***
              if (destination.amsl_flag)
              dest_alt=four_thirds_earth+destination.alt;
              else
             ***/
            dest_alt=four_thirds_earth+destination.alt+path->elevation[y];

            dest_alt2=dest_alt*dest_alt;
            source_alt2=source_alt*source_alt;

            /* Calculate the cosine of the elevation of
               the receiver as seen by the transmitter. */

            cos_xmtr_angle=((source_alt2)+(distance*distance)-(dest_alt2))/(2.0*source_alt*distance);

            if (got_elevation_pattern)
            {
                /* If an antenna elevation pattern is available, the
                   following code determines the elevation angle to
                   the first obstruction along the path-> */

                for (x=2, block=0; x<y && block==0; x++)
                {
                    distance=5280.0*(path->distance[y]-path->distance[x]);
                    test_alt=four_thirds_earth+path->elevation[x];

                    /* Calculate the cosine of the elevation
                       angle of the terrain (test point)
                       as seen by the transmitter. */

                    cos_test_angle=((source_alt2)+(distance*distance)-(test_alt*test_alt))/(2.0*source_alt*distance);

                    /* Compare these two angles to determine if
                       an obstruction exists.  Since we're comparing
                       the cosines of these angles rather than
                       the angles themselves, the sense of the
                       following "if" statement is reversed from
                       what it would be if the angles themselves
                       were compared. */

                    if (cos_xmtr_angle>=cos_test_angle)
                        block=1;
                }

                /* At this point, we have the elevation angle
                   to the first obstruction (if it exists). */
            }

            /* Determine path loss for each point along
               the path using ITWOM's point_to_point mode
               starting at x=2 (number_of_points = 1), the
               shortest distance terrain can play a role in
               path loss. */

            elev[0]=(elev_t)(y-1);	/* (number of points - 1) */

            /* Distance between elevation samples */

            elev[1]=(elev_t)(METERS_PER_MILE*(path->distance[y]-path->distance[y-1]));

            if (itwom)
                point_to_point(elev,source.alt*METERS_PER_FOOT, 
                               destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                               LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                               LR.radio_climate, LR.pol, LR.conf, LR.rel, &loss,
                               strmode, &errnum);
            else
                point_to_point_ITM(elev,source.alt*METERS_PER_FOOT, 
                                   destination.alt*METERS_PER_FOOT, LR.eps_dielect,
                                   LR.sgm_conductivity, LR.eno_ns_surfref, LR.frq_mhz,
                                   LR.radio_climate, LR.pol, LR.conf, LR.rel, &loss,
                                   strmode, &errnum);

            if (block)
                elevation=((acos(cos_test_angle))/DEG2RAD)-90.0;
            else
                elevation=((acos(cos_xmtr_angle))/DEG2RAD)-90.0;

            /* Integrate the antenna's radiation
               pattern into the overall path loss. */

            x=(int)rint(10.0*(10.0-elevation));

            if (x>=0 && x<=1000)
            {
                pattern=(double)LR.antenna_pattern[(int)azimuth][x];

                if (pattern > 0.0) {
                    patterndB=20.0*log10(pattern);
                } else if (pattern != NO_ANTENNA_DATA) {
                    loss-=9999;
                }
            }
            else
                patterndB=0.0;

            total_loss=loss-patterndB;

            if (metric)
                fprintf(fd,"%f\t%f\n",KM_PER_MILE*path->distance[y],total_loss);
            else
                fprintf(fd,"%f\t%f\n",path->distance[y],total_loss);

            if (total_loss>maxloss)
                maxloss=total_loss;

            if (total_loss<minloss)
                minloss=total_loss;
        }

        fclose(fd);

        free(elev);

        distance=Distance(source,destination);


        if (distance!=0.0)
        {
            free_space_loss=36.6+(20.0*log10(LR.frq_mhz))+(20.0*log10(distance));

            fprintf(fd2,"Free space path loss: %.2f dB\n",free_space_loss);
        }

        if (itwom)
            fprintf(fd2,"ITWOM Version %.1f path loss: %.2f dB\n",ITWOMVersion(),loss);
        else
            fprintf(fd2,"Longley-Rice path loss: %.2f dB\n",loss);

        if (free_space_loss!=0.0)
            fprintf(fd2,"Attenuation due to terrain shielding: %.2f dB\n",loss-free_space_loss);

        if (patterndB!=0.0)
            fprintf(fd2,"Total path loss including %s antenna pattern: %.2f dB\n",source.name,total_loss);

        if (LR.erp!=0.0)
        {
            field_strength=(139.4+(20.0*log10(LR.frq_mhz))-total_loss)+(10.0*log10(LR.erp/1000.0));

            /* dBm is referenced to EIRP */

            rxp=eirp/(pow(10.0,(total_loss/10.0)));
            dBm=10.0*(log10(rxp*1000.0));
            power_density=(eirp/(pow(10.0,(total_loss-free_space_loss)/10.0)));
            /* divide by 4*PI*distance_in_meters squared */
            power_density/=(4.0*PI*distance*distance*2589988.11);

            fprintf(fd2,"Field strength at %s: %.2f dBuV/meter\n", destination.name,field_strength);
            fprintf(fd2,"Signal power level at %s: %+.2f dBm\n",destination.name,dBm);
            fprintf(fd2,"Signal power density at %s: %+.2f dBW per square meter\n",destination.name,10.0*log10(power_density));
            voltage=1.0e6*sqrt(50.0*(eirp/(pow(10.0,(total_loss-2.14)/10.0))));
            fprintf(fd2,"Voltage across a 50 ohm dipole at %s: %.2f uV (%.2f dBuV)\n",destination.name,voltage,20.0*log10(voltage));

            voltage=1.0e6*sqrt(75.0*(eirp/(pow(10.0,(total_loss-2.14)/10.0))));
            fprintf(fd2,"Voltage across a 75 ohm dipole at %s: %.2f uV (%.2f dBuV)\n",destination.name,voltage,20.0*log10(voltage));
        }

        fprintf(fd2,"Mode of propagation: ");

        if (itwom)
        {
            if (strcmp(strmode,"L-o-S")==0)
                fprintf(fd2,"Line of Sight\n");

            if (strncmp(strmode,"1_Hrzn",6)==0)
                fprintf(fd2,"Single Horizon ");

            if (strncmp(strmode,"2_Hrzn",6)==0)
                fprintf(fd2,"Double Horizon ");

            y= (int)strlen(strmode);

            if (y>19)
                y=19;

            for (x=6; x<y; x++)
                propbuf[x-6]=strmode[x];

            propbuf[x]=0;

            if (strncmp(propbuf,"_Diff",5)==0)
                fprintf(fd2,"Diffraction Dominant\n");

            if (strncmp(propbuf,"_Tropo",6)==0)
                fprintf(fd2,"Troposcatter Dominant\n");

            if (strncmp(propbuf,"_Peak",5)==0)
                fprintf(fd2,"RX at Peak Terrain Along Path\n");

            fprintf(fd2,"ITWOM error number: %d",errnum);
        }
        else
        {
            fprintf(fd2,"%s\n",strmode);
            fprintf(fd2,"Longley-Rice model error number: %d",errnum);
        }


        switch (errnum)
        {
            case 0:
                fprintf(fd2," (No error)\n");
                break;

            case 1:
                fprintf(fd2,"\n  Warning: Some parameters are nearly out of range.\n");
                fprintf(fd2,"  Results should be used with caution.\n");
                break;

            case 2:
                fprintf(fd2,"\n  Note: Default parameters have been substituted for impossible ones.\n");
                break;

            case 3:
                fprintf(fd2,"\n  Warning: A combination of parameters is out of range.\n");
                fprintf(fd2,"  Results are probably invalid.\n");
                break;

            default:
                fprintf(fd2,"\n  Warning: Some parameters are out of range.\n");
                fprintf(fd2,"  Results are probably invalid.\n");
        }

        fprintf(fd2,"\n%s\n\n",dashes);
    }

    fprintf(stdout,"\nPath Loss Report written to: \"%s\"\n",report_name);
    fflush(stdout);

    ObstructionAnalysis(source, destination, LR.frq_mhz, fd2);

    fclose(fd2);

    /* Skip plotting the graph if ONLY a path-loss report is needed. */

    if (graph_it)
    {
        if (name[0]=='.')
        {
            /* Default filename and output file type */
            strcpy(basename,"profile");
            strcpy(ext,"png");
        }
        else
        {
            /* Extract extension and terminal type from "name" */
            strncpy(basename, name, MAX_PATH_LEN-1);
            stripExtension(basename, ext, 20);
            convertBackslashes(basename);
            if (strlen(ext) == 0)
            {
                strcpy(ext,"png");
            }
        }
        strcpy(term,ext);

        /* Either .ps or .postscript may be used
           as an extension for postscript output. */

        if (strncmp(term,"postscript",10)==0)
            strcpy(ext,"ps0");

        else if (strncmp(ext,"ps",2)==0)
            strcpy(term,"postscript enhanced color0");

        fd=fopen("splat.gp","w");

        fprintf(fd,"set grid\n");
        fprintf(fd,"set yrange [%2.3f to %2.3f]\n", minloss, maxloss);
        fprintf(fd,"set encoding iso_8859_1\n");
        fprintf(fd,"set term %s\n",term);
        fprintf(fd,"set title \"%s Loss Profile Along Path Between %s and %s (%.2f%c azimuth)\"\n",SPLAT_NAME, destination.name, source.name, Azimuth(destination,source),176);

        if (metric)
            fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f kilometers)\"\n",destination.name,source.name,KM_PER_MILE*Distance(destination,source));
        else
            fprintf(fd,"set xlabel \"Distance Between %s and %s (%.2f miles)\"\n",destination.name,source.name,Distance(destination,source));

        if (got_azimuth_pattern || got_elevation_pattern)
            fprintf(fd,"set ylabel \"Total Path Loss (including TX antenna pattern) (dB)");
        else
        {
            if (itwom)
                fprintf(fd,"set ylabel \"ITWOM Version %.1f Path Loss (dB)",ITWOMVersion());
            else
                fprintf(fd,"set ylabel \"Longley-Rice Path Loss (dB)");
        }

        fprintf(fd,"\"\nset output \"%s.%s\"\n",basename,ext);
        fprintf(fd,"plot \"profile.gp\" title \"Path Loss\" with lines\n");

        fclose(fd);

        x=system("gnuplot splat.gp");

        if (x!=-1)
        {
            if (!gpsav)
            {
                unlink("splat.gp");
                unlink("profile.gp");
                unlink("reference.gp");
            }

            fprintf(stdout,"Path loss plot written to: \"%s.%s\"\n",basename,ext);
            fflush(stdout);
        }

        else
            fprintf(stderr,"\n*** ERROR: Error occurred invoking gnuplot!\n");
    }

    if (x!=-1 && !gpsav)
        unlink("profile.gp");

    DestroyPath(path);

}

/* Generates a site analysis report.
 */
void SiteReport(Site xmtr)
{
    char	report_name[MAX_PATH_LEN];
    double	terrain;
    int	x, azi;
    FILE	*fd;

    snprintf(report_name,MAX_PATH_LEN,"%s-site_report.txt",xmtr.name);

    for (x=0; report_name[x]!=0; x++)
        if (report_name[x]==32 || report_name[x]==17 || report_name[x]==92 || report_name[x]==42 || report_name[x]==47)
            report_name[x]='_';

    fd=fopen(report_name,"w");

    fprintf(fd,"\n\t--==[ %s v%s Site Analysis Report For: %s ]==--\n\n",SPLAT_NAME, SPLAT_VERSION, xmtr.name);

    fprintf(fd,"%s\n\n",dashes);

    if (xmtr.lat>=0.0)
    {
        fprintf(fd,"Site location: %.4f North / %.4f West",xmtr.lat, xmtr.lon);
        fprintf(fd, " (%s N / ",dec2dms(xmtr.lat));
    }

    else
    {
        fprintf(fd,"Site location: %.4f South / %.4f West",-xmtr.lat, xmtr.lon);
        fprintf(fd, " (%s S / ",dec2dms(xmtr.lat));
    }

    fprintf(fd, "%s W)\n",dec2dms(xmtr.lon));

    if (metric)
    {
        if (xmtr.amsl_flag)
            fprintf(fd,"Antenna height: %.2f meters AMSL\n",METERS_PER_FOOT*xmtr.alt);
        else
        {
            fprintf(fd,"Ground elevation: %.2f meters AMSL\n",METERS_PER_FOOT*GetElevation(xmtr));
            fprintf(fd,"Antenna height: %.2f meters AGL / %.2f meters AMSL\n",METERS_PER_FOOT*xmtr.alt, METERS_PER_FOOT*(xmtr.alt+GetElevation(xmtr)));
        }
    }

    else
    {
        if (xmtr.amsl_flag)
            fprintf(fd,"Antenna height: %.2f feet AMSL\n",xmtr.alt);
        else
        {
            fprintf(fd,"Ground elevation: %.2f feet AMSL\n",GetElevation(xmtr));
            fprintf(fd,"Antenna height: %.2f feet AGL / %.2f feet AMSL\n",xmtr.alt, xmtr.alt+GetElevation(xmtr));
        }
    }

    terrain=haat(xmtr);

    if (terrain>-4999.0)
    {
        if (metric)
            fprintf(fd,"Antenna height above average terrain: %.2f meters\n\n",METERS_PER_FOOT*terrain);
        else
            fprintf(fd,"Antenna height above average terrain: %.2f feet\n\n",terrain);

        /* Display the average terrain between 2 and 10 miles
           from the transmitter site at azimuths of 0, 45, 90,
           135, 180, 225, 270, and 315 degrees. */

        for (azi=0; azi<=315; azi+=45)
        {
            fprintf(fd,"Average terrain at %3d degrees azimuth: ",azi);
            terrain=AverageTerrain(xmtr,(double)azi,2.0,10.0);

            if (terrain>-4999.0)
            {
                if (metric)
                    fprintf(fd,"%.2f meters AMSL\n",METERS_PER_FOOT*terrain);
                else
                    fprintf(fd,"%.2f feet AMSL\n",terrain);
            }

            else
                fprintf(fd,"No terrain\n");
        }
    }

    fprintf(fd,"\n%s\n\n",dashes);
    fclose(fd);
    fprintf(stdout,"\nSite analysis report written to: \"%s\"\n",report_name);
}

/* Loads the SDF files required to cover the limits of the region
 * specified.
 */ 
void LoadTopoData(int max_lon, int min_lon, int max_lat, int min_lat)
{
    int x, y, width, ymin, ymax;

    width=ReduceAngle(max_lon-min_lon);

    if ((max_lon-min_lon)<=180.0)
    {
        for (y=0; y<=width; y++)
            for (x=min_lat; x<=max_lat; x++)
            {
                ymin=(int)(min_lon+(double)y);

                while (ymin<0)
                    ymin+=360;

                while (ymin>=360)
                    ymin-=360;

                ymax=ymin+1;

                while (ymax<0)
                    ymax+=360;

                while (ymax>=360)
                    ymax-=360;

                LoadSDF(x, x+1, ymin, ymax, (ippd==3600));
            }
    }

    else
    {
        for (y=0; y<=width; y++)
            for (x=min_lat; x<=max_lat; x++)
            {
                ymin=max_lon+y;

                while (ymin<0)
                    ymin+=360;

                while (ymin>=360)
                    ymin-=360;

                ymax=ymin+1;

                while (ymax<0)
                    ymax+=360;

                while (ymax>=360)
                    ymax-=360;

                LoadSDF(x, x+1, ymin, ymax, (ippd==3600));
            }
    }
}

/* Reads a SPLAT! alphanumeric output file (-ani option) for analysis
 * and/or map generation.
 */
int LoadANO(char *filename)
{
    int	error=0, max_west, min_west, max_north, min_north;
    char	buf[MAX_LINE_LEN], *pointer=NULL;
    double	latitude=0.0, longitude=0.0, azimuth=0.0, elevation=0.0,
            ano=0.0;
    FILE	*fd;

    fd=fopen(filename,"r");

    if (fd!=NULL)
    {
        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(buf,"%d, %d",&max_west, &min_west);

        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        sscanf(buf,"%d, %d",&max_north, &min_north);

        fgets(buf,MAX_LINE_LEN,fd);
        pointer=strchr(buf,';');

        if (pointer!=NULL)
            *pointer=0;

        LoadTopoData(max_west-1, min_west, max_north-1, min_north);

        fprintf(stdout,"\nReading \"%s\"... ",filename);
        fflush(stdout);

        fgets(buf,MAX_LINE_LEN,fd);
        sscanf(buf,"%lf, %lf, %lf, %lf, %lf",&latitude, &longitude, &azimuth, &elevation, &ano);

        while (feof(fd)==0)
        {
            if (LR.erp==0.0)
            {
                /* Path loss */

                if (contour_threshold==0 || (fabs(ano)<=(double)contour_threshold))
                {
                    ano=fabs(ano);

                    if (ano>255.0)
                        ano=255.0;

                    PutSignal(latitude,longitude,((unsigned char)round(ano)));
                }
            }

            if (LR.erp!=0.0 && dbm)
            {
                /* signal power level in dBm */

                if (contour_threshold==0 || (ano>=(double)contour_threshold))
                {
                    ano=200.0+rint(ano);

                    if (ano<0.0)
                        ano=0.0;

                    if (ano>255.0)
                        ano=255.0;

                    PutSignal(latitude,longitude,((unsigned char)round(ano)));
                }
            }

            if (LR.erp!=0.0 && !dbm)
            {
                /* field strength dBuV/m */

                if (contour_threshold==0 || (ano>=(double)contour_threshold))
                {
                    ano=100.0+rint(ano);

                    if (ano<0.0)
                        ano=0.0;

                    if (ano>255.0)
                        ano=255.0;

                    PutSignal(latitude,longitude,((unsigned char)round(ano)));
                }
            }

            fgets(buf,MAX_LINE_LEN,fd);
            sscanf(buf,"%lf, %lf, %lf, %lf, %lf",&latitude, &longitude, &azimuth, &elevation, &ano);
        }

        fclose(fd);

        fprintf(stdout," Done!\n");
        fflush(stdout);
    }

    else
        error=1;

    return error;
}

/* Writes a KML file plotting the path between source and destination.
 */
void WriteKML(Site source, Site destination)
{
    int	x, y;
    char	block;
    char report_name[MAX_PATH_LEN-4], kmlfile[MAX_PATH_LEN];
#if DO_KMZ
    char kmzfile[MAX_PATH_LEN];
#endif
    double	distance, rx_alt, tx_alt, cos_xmtr_angle,
            azimuth, cos_test_angle, test_alt;
    FILE	*fd=NULL;

    Path *path = AllocatePath();
    if (!path) {
        fprintf(stderr,"\n*** ERROR: Couldn't allocate memory for Path");
        return;
    }
    ReadPath(source,destination,path);

    snprintf(report_name,MAX_PATH_LEN-4,"%s-to-%s",source.name,destination.name);

    for (x=0; report_name[x]!=0; x++)
        if (report_name[x]==32 || report_name[x]==17 || report_name[x]==92 || report_name[x]==42 || report_name[x]==47)
            report_name[x]='_';

    snprintf(kmlfile, MAX_PATH_LEN, "%s.kml", report_name);
    fd=fopen(report_name,"w");

    fprintf(fd,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fd,"<kml xmlns=\"http://earth.google.com/kml/2.0\">\n");
    fprintf(fd,"<!-- Generated by %s Version %s -->\n",SPLAT_NAME, SPLAT_VERSION);
    fprintf(fd,"<Folder>\n");
    fprintf(fd,"<name>SPLAT! Path</name>\n");
    fprintf(fd,"<open>1</open>\n");
    fprintf(fd,"<description>Path Between %s and %s</description>\n",source.name,destination.name);

    fprintf(fd,"<Placemark>\n");
    fprintf(fd,"	<name>%s</name>\n",source.name);
    fprintf(fd,"	<description>\n");
    fprintf(fd,"	   Transmit Site\n");

    if (source.lat>=0.0)
        fprintf(fd,"	   <BR>%s North</BR>\n",dec2dms(source.lat));
    else
        fprintf(fd,"	   <BR>%s South</BR>\n",dec2dms(source.lat));

    fprintf(fd,"	   <BR>%s West</BR>\n",dec2dms(source.lon));

    azimuth=Azimuth(source,destination);
    distance=Distance(source,destination);

    if (metric)
        fprintf(fd,"	   <BR>%.2f km",distance*KM_PER_MILE);
    else
        fprintf(fd,"	   <BR>%.2f miles",distance);

    fprintf(fd," to %s</BR>\n	   <BR>toward an azimuth of %.2f%c</BR>\n",destination.name,azimuth,176);

    fprintf(fd,"	</description>\n");
    fprintf(fd,"	<visibility>1</visibility>\n");
    fprintf(fd,"	<Style>\n");
    fprintf(fd,"	  <IconStyle>\n");
    fprintf(fd,"		<Icon>\n");
    fprintf(fd,"		  <href>root://icons/palette-5.png</href>\n");
    fprintf(fd,"		  <x>224</x>\n");
    fprintf(fd,"		  <y>224</y>\n");
    fprintf(fd,"		  <w>32</w>\n");
    fprintf(fd,"		  <h>32</h>\n");
    fprintf(fd,"		</Icon>\n");
    fprintf(fd,"	  </IconStyle>\n");
    fprintf(fd,"	</Style>\n");
    fprintf(fd,"	<Point>\n");
    fprintf(fd,"	  <extrude>1</extrude>\n");
    fprintf(fd,"	  <altitudeMode>relativeToGround</altitudeMode>\n");
    fprintf(fd,"	  <coordinates>%f,%f,30</coordinates>\n",(source.lon<180.0?-source.lon:360.0-source.lon),source.lat);
    fprintf(fd,"	</Point>\n");
    fprintf(fd,"</Placemark>\n");

    fprintf(fd,"<Placemark>\n");
    fprintf(fd,"	<name>%s</name>\n",destination.name);
    fprintf(fd,"	<description>\n");
    fprintf(fd,"	   Receive Site\n");

    if (destination.lat>=0.0)
        fprintf(fd,"	   <BR>%s North</BR>\n",dec2dms(destination.lat));
    else
        fprintf(fd,"	   <BR>%s South</BR>\n",dec2dms(destination.lat));

    fprintf(fd,"	   <BR>%s West</BR>\n",dec2dms(destination.lon));

    if (metric)
        fprintf(fd,"	   <BR>%.2f km",distance*KM_PER_MILE);
    else
        fprintf(fd,"	   <BR>%.2f miles",distance);

    fprintf(fd," to %s</BR>\n	   <BR>toward an azimuth of %.2f%c</BR>\n",source.name,Azimuth(destination,source),176);

    fprintf(fd,"	</description>\n");
    fprintf(fd,"	<visibility>1</visibility>\n");
    fprintf(fd,"	<Style>\n");
    fprintf(fd,"	  <IconStyle>\n");
    fprintf(fd,"		<Icon>\n");
    fprintf(fd,"		  <href>root://icons/palette-5.png</href>\n");
    fprintf(fd,"		  <x>224</x>\n");
    fprintf(fd,"		  <y>224</y>\n");
    fprintf(fd,"		  <w>32</w>\n");
    fprintf(fd,"		  <h>32</h>\n");
    fprintf(fd,"		</Icon>\n");
    fprintf(fd,"	  </IconStyle>\n");
    fprintf(fd,"	</Style>\n");
    fprintf(fd,"	<Point>\n");
    fprintf(fd,"	  <extrude>1</extrude>\n");
    fprintf(fd,"	  <altitudeMode>relativeToGround</altitudeMode>\n");
    fprintf(fd,"	  <coordinates>%f,%f,30</coordinates>\n",(destination.lon<180.0?-destination.lon:360.0-destination.lon),destination.lat);
    fprintf(fd,"	</Point>\n");
    fprintf(fd,"</Placemark>\n");

    fprintf(fd,"<Placemark>\n");
    fprintf(fd,"<name>Point-to-Point Path</name>\n");
    fprintf(fd,"  <visibility>1</visibility>\n");
    fprintf(fd,"  <open>0</open>\n");
    fprintf(fd,"  <Style>\n");
    fprintf(fd,"	<LineStyle>\n");
    fprintf(fd,"	  <color>7fffffff</color>\n");
    fprintf(fd,"	</LineStyle>\n");
    fprintf(fd,"	<PolyStyle>\n");
    fprintf(fd,"	   <color>7fffffff</color>\n");
    fprintf(fd,"	</PolyStyle>\n");
    fprintf(fd,"  </Style>\n");
    fprintf(fd,"  <LineString>\n");
    fprintf(fd,"	<extrude>1</extrude>\n");
    fprintf(fd,"	<tessellate>1</tessellate>\n");
    fprintf(fd,"	<altitudeMode>relativeToGround</altitudeMode>\n");
    fprintf(fd,"	<coordinates>\n");

    for (x=0; x<path->length; x++)
        fprintf(fd,"	  %f,%f,5\n",(path->lon[x]<180.0?-path->lon[x]:360.0-path->lon[x]),path->lat[x]);

    fprintf(fd,"	</coordinates>\n");
    fprintf(fd,"   </LineString>\n");
    fprintf(fd,"</Placemark>\n");

    fprintf(fd,"<Placemark>\n");
    fprintf(fd,"<name>Line-of-Sight Path</name>\n");
    fprintf(fd,"  <visibility>1</visibility>\n");
    fprintf(fd,"  <open>0</open>\n");
    fprintf(fd,"  <Style>\n");
    fprintf(fd,"	<LineStyle>\n");
    fprintf(fd,"	  <color>ff00ff00</color>\n");
    fprintf(fd,"	</LineStyle>\n");
    fprintf(fd,"	<PolyStyle>\n");
    fprintf(fd,"	   <color>7f00ff00</color>\n");
    fprintf(fd,"	</PolyStyle>\n");
    fprintf(fd,"  </Style>\n");
    fprintf(fd,"  <LineString>\n");
    fprintf(fd,"	<extrude>1</extrude>\n");
    fprintf(fd,"	<tessellate>1</tessellate>\n");
    fprintf(fd,"	<altitudeMode>relativeToGround</altitudeMode>\n");
    fprintf(fd,"	<coordinates>\n");

    /* Walk across the "path", indentifying obstructions along the way */

    for (y=0; y<path->length; y++)
    {
        distance=5280.0*path->distance[y];

        if (source.amsl_flag)
            tx_alt=earthradius+source.alt;
        else
            tx_alt=earthradius+source.alt+path->elevation[0];

        /***
          if (destination.amsl_flag)
          rx_alt=earthradius+destination.alt;
          else
         ***/
        rx_alt=earthradius+destination.alt+path->elevation[y];

        /* Calculate the cosine of the elevation of the
           transmitter as seen at the temp rx point. */

        cos_xmtr_angle=((rx_alt*rx_alt)+(distance*distance)-(tx_alt*tx_alt))/(2.0*rx_alt*distance);

        for (x=y, block=0; x>=0 && block==0; x--)
        {
            distance=5280.0*(path->distance[y]-path->distance[x]);
            test_alt=earthradius+path->elevation[x];

            cos_test_angle=((rx_alt*rx_alt)+(distance*distance)-(test_alt*test_alt))/(2.0*rx_alt*distance);

            /* Compare these two angles to determine if
               an obstruction exists.  Since we're comparing
               the cosines of these angles rather than
               the angles themselves, the following "if"
               statement is reversed from what it would
               be if the actual angles were compared. */

            if (cos_xmtr_angle>=cos_test_angle)
                block=1;
        }

        if (block)
            fprintf(fd,"	  %f,%f,-30\n",(path->lon[y]<180.0?-path->lon[y]:360.0-path->lon[y]),path->lat[y]);
        else
            fprintf(fd,"	  %f,%f,5\n",(path->lon[y]<180.0?-path->lon[y]:360.0-path->lon[y]),path->lat[y]);
    }

    fprintf(fd,"	</coordinates>\n");
    fprintf(fd,"  </LineString>\n");
    fprintf(fd,"</Placemark>\n");

    fprintf(fd,"	<LookAt>\n");
    fprintf(fd,"	  <longitude>%f</longitude>\n",(source.lon<180.0?-source.lon:360.0-source.lon));
    fprintf(fd,"	  <latitude>%f</latitude>\n",source.lat);
    fprintf(fd,"	  <range>300.0</range>\n");
    fprintf(fd,"	  <tilt>45.0</tilt>\n");
    fprintf(fd,"	  <heading>%f</heading>\n",azimuth);
    fprintf(fd,"	</LookAt>\n");

    fprintf(fd,"</Folder>\n");
    fprintf(fd,"</kml>\n");

    fclose(fd);


#if DO_KMZ
    snprintf(kmzfile, MAX_PATH_LEN, "%s.kmz", report_name);
    bool success = false;

    struct zip_t *zip = zip_open(kmzfile, ZIP_DEFAULT_COMPRESSION_LEVEL, 'w');
    if (zip) {
        if (zip_entry_open(zip, kmlfile) == 0) {
            if (zip_entry_fwrite(zip, kmlfile) == 0) {
                success = true;
            }
            zip_entry_close(zip);
        }
        zip_close(zip);
    }

    if (success) {
        unlink(kmlfile);
        fprintf(stdout, "\nKMZ file written to: \"%s\"\n", kmzfile);
    } else {
        unlink(kmzfile);
        fprintf(stdout, "\nCouldn't create KMZ file. KML file written to: \"%s\"\n", kmlfile);
    }
#else
    fprintf(stdout, "\nKML file written to: \"%s\"", kmlfile);
#endif

    fflush(stdout);

    DestroyPath(path);
}

int main(int argc, char *argv[])
{
    int     errcount = 0;
    int		min_lat, min_lon, max_lat, max_lon,
            rxlat, rxlon, txlat, txlon, west_min, west_max,
            north_min, north_max;
    unsigned int txsites=0, max_txsites=0, cities=0, bfs=0;
    bool	coverage=false, LRmap=false, terrain_plot=false,
            elevation_plot=false, height_plot=false, norm_height_plot=false,
            map=false, longley_plot=false;
    bool    topomap=false, geo=false, kml=false, pt2pt_mode=false,
            area_mode=false, ngs=false, nolospath=false,
            nositereports=false, fresnel_plot=true, command_line_log=false;

#ifdef HAVE_LIBPNG
    ImageType imagetype=IMAGETYPE_PNG;
#else
    ImageType imagetype=IMAGETYPE_PPM;
#endif
    unsigned char imagetype_set = 0;

    char    header[80];
    char	mapfile[MAX_PATH_LEN], city_file[5][MAX_PATH_LEN];
    char	longley_file[MAX_PATH_LEN-8], terrain_file[MAX_PATH_LEN-8],
            elevation_file[MAX_PATH_LEN-8], height_file[MAX_PATH_LEN-8],
            norm_height_file[MAX_PATH_LEN-8];
    char	buf[MAX_PATH_LEN], rxfile[MAX_PATH_LEN], 
            txfile[MAX_PATH_LEN], boundary_file[5][MAX_PATH_LEN],
            udt_file[MAX_PATH_LEN], ani_filename[MAX_PATH_LEN],
            ano_filename[MAX_PATH_LEN], logfile[MAX_PATH_LEN];

    char    *env = NULL;
    char    rxsite = 0;

    double		altitude=0.0, altitudeLR=0.0, tx_range=0.0,
                rx_range=0.0, deg_range=0.0, deg_limit=0.0,
                deg_range_lon, er_mult;

    bool	multithread = true;

    Site	tx_site[32], rx_site;

    FILE		*fd;

    unsigned long long systemMemory = getTotalSystemMemory();

    strncpy(dashes,"---------------------------------------------------------------------------\0",76);

    /* backwards compatibility - set to HD mode if invoked as splat-hd */
    if (strstr(argv[0], "-hd")!=NULL) {
        appmode = APPMODE_HD;
    }

    if (argc==1)
    {
        fprintf(stdout,"\n\t\t --==[ %s%s v%s Available Options... ]==--\n\n",
                SPLAT_NAME, (appmode==APPMODE_HD?" HD":""), SPLAT_VERSION);

        fprintf(stdout,"       -t txsite(s).qth (max of 4 with -c, max of 30 with -L)\n");
        fprintf(stdout,"       -r rxsite.qth\n");
        fprintf(stdout,"       -c plot LOS coverage of TX(s) with an RX antenna at X feet/meters AGL\n");
        fprintf(stdout,"       -L plot path loss map of TX based on an RX at X feet/meters AGL\n");
        fprintf(stdout,"       -s filename(s) of city/site file(s) to import (5 max)\n");
        fprintf(stdout,"       -b filename(s) of cartographic boundary file(s) to import (5 max)\n");
        fprintf(stdout,"       -p filename of terrain profile graph to plot\n");
        fprintf(stdout,"       -e filename of terrain elevation graph to plot\n");
        fprintf(stdout,"       -h filename of terrain height graph to plot\n");
        fprintf(stdout,"       -H filename of normalized terrain height graph to plot\n");
        fprintf(stdout,"       -l filename of path loss graph to plot\n");
        fprintf(stdout,"       -o filename of topographic map to generate (without suffix)\n");
        fprintf(stdout,"       -u filename of user-defined terrain file to import\n");
        fprintf(stdout,"       -d sdf file directory path (overrides path in ~/.splat_path file)\n");
        fprintf(stdout,"       -m earth radius multiplier\n");
        fprintf(stdout,"       -n do not plot LOS paths in output maps\n");
        fprintf(stdout,"       -N do not produce unnecessary site or obstruction reports\n");
        fprintf(stdout,"       -f frequency for Fresnel zone calculation (MHz)\n");
        fprintf(stdout,"       -R modify default range for -c or -L (miles/kilometers)\n");
        fprintf(stdout,"       -v N verbosity level. Default is 1. Set to 0 to quiet everything.\n");
        fprintf(stdout,"      -st use a single CPU thread (classic mode)\n");
        fprintf(stdout,"      -hd Use High Definition mode. Requires 1-deg SDF files.\n");
        fprintf(stdout,"      -sc display smooth rather than quantized contour levels\n");
        fprintf(stdout,"      -db threshold beyond which contours will not be displayed\n");
        fprintf(stdout,"      -nf do not plot Fresnel zones in height plots\n");
        fprintf(stdout,"      -fz Fresnel zone clearance percentage (default = 60)\n");
        fprintf(stdout,"      -gc ground clutter height (feet/meters)\n");
        fprintf(stdout,"     -jpg when generating maps, create jpgs instead of pngs or ppms\n");
#ifdef HAVE_LIBPNG
        fprintf(stdout,"     -ppm when generating maps, create ppms instead of pngs or jpgs\n");
#endif
        fprintf(stdout,"     -ngs display greyscale topography as white in .ppm files\n"); 
        fprintf(stdout,"     -erp override ERP in .lrp file (Watts)\n");
        fprintf(stdout,"     -ano name of alphanumeric output file\n");
        fprintf(stdout,"     -ani name of alphanumeric input file\n");
#ifndef _WIN32
        fprintf(stdout,"     -udt name of user defined terrain input file\n");
#endif
        fprintf(stdout,"     -kml generate Google Earth (.kml) compatible output\n");
        fprintf(stdout,"     -geo generate an Xastir .geo georeference file (with .ppm output)\n");
        fprintf(stdout,"     -dbm plot signal power level contours rather than field strength\n");
        fprintf(stdout,"     -log copy command line string to this output file\n");
        fprintf(stdout,"   -gpsav preserve gnuplot temporary working files after SPLAT! execution\n");
        fprintf(stdout,"   -itwom invoke the ITWOM model instead of using Longley-Rice\n");
        fprintf(stdout,"  -metric employ metric rather than imperial units for all user I/O\n");

        fprintf(stdout,"See the documentation for more details.\n\n");

        fprintf(stdout, "Using ITWOM Version %.1f.\n\n",ITWOMVersion());

        fprintf(stdout,"This compilation of %s supports analysis over a region of %dx%d square degree%s.\n",
                SPLAT_NAME, maxpagesides, maxpagesides, (maxpagesides>1)?"s":"");

        unsigned long pathmem;
        unsigned long demmem;
        pathmem = (WorkQueue::maxWorkers()*SizeofPath(APPMODE_NORMAL, maxpagesides))/1024;
        demmem = (maxpagesides*maxpagesides*SizeofDEM(APPMODE_NORMAL))/1024;
        fprintf(stdout,"At the %dx%d maximum in normal mode it will require %lu MB of RAM.\n",
                maxpagesides, maxpagesides, (pathmem+demmem)/1024);
        fprintf(stdout,"Each Path can be up to %lu KB.\n\n", pathmem);
        /*
           fprintf(stdout,"Each DEM: %lu KB\n", demmem);
           fprintf(stdout,"\n");
         */

        pathmem = (WorkQueue::maxWorkers()*SizeofPath(APPMODE_HD, maxpagesides))/1024;
        demmem = (maxpagesides*maxpagesides*SizeofDEM(APPMODE_HD))/1024;
        fprintf(stdout,"At the %dx%d maximum in HD mode it will require %lu MB of RAM.\n",
                maxpagesides, maxpagesides, (pathmem+demmem)/1024);
        fprintf(stdout,"Each Path can be up to %lu KB.\n\n", pathmem);
        /*
           fprintf(stdout,"Each DEM: %lu KB\n", demmem);
           fprintf(stdout,"\n");
         */

        fprintf(stdout,"You have %lu MB of RAM available.\n", (unsigned long)(systemMemory/(1024*1024)) );
        fprintf(stdout,"You have %d cores/threads available.\n", WorkQueue::maxWorkers());

        fflush(stdout);

        return 0;
    }

    kml=0;
    geo=0;
    metric=0;
    rxfile[0]=0;
    txfile[0]=0;
    buf[0]=0;
    mapfile[0]=0;
    clutter=0.0;
    forced_erp=-1.0;
    forced_freq=0.0;
    elevation_file[0]=0;
    terrain_file[0]=0;
    sdf_path[0]=0;
    home_sdf_path[0]=0;
    udt_file[0]=0;
    max_txsites=30;
    fzone_clearance=0.6;
    contour_threshold=0;
    rx_site.lat=91.0;
    rx_site.lon=361.0;
    longley_file[0]=0;
    ano_filename[0]=0;
    ani_filename[0]=0;
    earthradius=EARTHRADIUS;

    for (int x=0; x<4; x++)
    {
        tx_site[x].lat=91.0;
        tx_site[x].lon=361.0;
    }

    /* Scan for command line arguments */

    for (int x=1, z=0; x<argc; x++)
    {
        if (strcmp(argv[x],"-R")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&max_range);

                if (max_range<0.0)
                    max_range=0.0;

                if (max_range>1000.0)
                    max_range=1000.0;
            }			 
        }

        if (strcmp(argv[x],"-m")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&er_mult);

                if (er_mult<0.1)
                    er_mult=1.0;

                if (er_mult>1.0e6)
                    er_mult=1.0e6;

                earthradius*=er_mult;
            }			 
        }

        if (strcmp(argv[x],"-v")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%d",&verbose);
            }			 
        }

        if (strcmp(argv[x],"-gc")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&clutter);

                if (clutter<0.0)
                    clutter=0.0;
            }			 
        }

        if (strcmp(argv[x],"-fz")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&fzone_clearance);

                if (fzone_clearance<0.0 || fzone_clearance>100.0)
                    fzone_clearance=60.0;

                fzone_clearance/=100.0;
            }
        }

        if (strcmp(argv[x],"-o")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
                strncpy(mapfile,argv[z],MAX_PATH_LEN-1);
            map=1;
        }

        if (strcmp(argv[x],"-log")==0)
        {
            z=x+1;

            logfile[0]=0;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
                strncpy(logfile,argv[z],MAX_PATH_LEN-1);

            command_line_log=1;
        }

#ifndef _WIN32
        if (strcmp(argv[x],"-udt")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
                strncpy(udt_file,argv[z],MAX_PATH_LEN-1);
        }
#endif

        if (strcmp(argv[x],"-c")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&altitude);
                map=1;
                coverage=1;
                area_mode=1;
                max_txsites=4;
            }
        }

        if (strcmp(argv[x],"-db")==0 || strcmp(argv[x],"-dB")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0]) /* A minus argument is legal here */
                sscanf(argv[z],"%d",&contour_threshold);
        }

        if (strcmp(argv[x],"-p")==0)
        { 
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                strncpy(terrain_file,argv[z],MAX_PATH_LEN-10);
                terrain_plot=1;
                pt2pt_mode=1;
            }
        }

        if (strcmp(argv[x],"-e")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                strncpy(elevation_file,argv[z],MAX_PATH_LEN-10);
                elevation_plot=1;
                pt2pt_mode=1;
            }
        }

        if (strcmp(argv[x],"-h")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                strncpy(height_file,argv[z],MAX_PATH_LEN-10);
                height_plot=1;
                pt2pt_mode=1;
            }
        }

        if (strcmp(argv[x],"-H")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                strncpy(norm_height_file,argv[z],MAX_PATH_LEN-10);
                norm_height_plot=1;
                pt2pt_mode=1;
            }
        }

#ifdef HAVE_LIBPNG
        if (strcmp(argv[x], "-ppm") == 0) {
            if (imagetype_set && imagetype != IMAGETYPE_PPM) {
                fprintf(stdout, "-jpg and -ppm are exclusive options, ignoring -ppm.\n");
            }
            else {
                imagetype = IMAGETYPE_PPM;
                imagetype_set = 1;
            }
        }
#endif
#ifdef HAVE_LIBJPEG
        if (strcmp(argv[x],"-jpg")==0) {
            if (imagetype_set && imagetype != IMAGETYPE_JPG) {
#ifdef HAVE_LIBPNG
                fprintf(stdout,"-jpg and -ppm are exclusive options, ignoring -jpg.\n");
#else
                fprintf(stdout,"-jpg and -png are exclusive options, ignoring -jpg.\n");
#endif
            } else {
                imagetype=IMAGETYPE_JPG;
                imagetype_set = 1;
            }
        }

#endif

        if (strcmp(argv[x],"-metric")==0)
            metric=1;

        if (strcmp(argv[x],"-gpsav")==0)
            gpsav=true;

        if (strcmp(argv[x],"-geo")==0)
            geo=1;

        if (strcmp(argv[x],"-kml")==0)
            kml=1;

        if (strcmp(argv[x],"-nf")==0)
            fresnel_plot=0;

        if (strcmp(argv[x],"-ngs")==0)
            ngs=1;

        if (strcmp(argv[x],"-n")==0)
            nolospath=1;

        if (strcmp(argv[x],"-dbm")==0)
            dbm=true;

        if (strcmp(argv[x],"-sc")==0)
            smooth_contours=true;

        if (strcmp(argv[x],"-st")==0)
            multithread=false;

        if (strcmp(argv[x],"-hd")==0) {
            appmode = APPMODE_HD;
        }

        if (strcmp(argv[x],"-itwom")==0)
            itwom=true;

        if (strcmp(argv[x],"-N")==0)
        {
            nolospath=1;
            nositereports=1;
        }

        if (strcmp(argv[x],"-d")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
                strncpy(sdf_path,argv[z],MAX_PATH_LEN-2);

            /* Ensure it has a trailing slash */
            size_t len = strlen(sdf_path);
            if (len > 0 && sdf_path[len-1]!=PATHSEP) {
                sdf_path[len] = PATHSEP;
                sdf_path[len + 1] = '\0';
            }
        }

        if (strcmp(argv[x],"-t")==0)
        {
            /* Read Transmitter Location */

            z=x+1;

            while (z<argc && argv[z][0] && argv[z][0]!='-' && txsites<30)
            {
                strncpy(txfile,argv[z],MAX_PATH_LEN-1);
                tx_site[txsites]=LoadQTH(txfile);
                txsites++;
                z++;
            }

            z--;
        }

        if (strcmp(argv[x],"-L")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&altitudeLR);
                map=1;
                LRmap=1;
                area_mode=1;

                if (coverage)
                    fprintf(stdout,"c and L are exclusive options, ignoring L.\n");
            }
        }

        if (strcmp(argv[x],"-l")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                strncpy(longley_file,argv[z],MAX_PATH_LEN-10);
                longley_plot=1;
                pt2pt_mode=1;
            }
        }

        if (strcmp(argv[x],"-r")==0)
        {
            /* Read Receiver Location */

            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                strncpy(rxfile,argv[z],MAX_PATH_LEN-1);
                rx_site=LoadQTH(rxfile);
                rxsite=1;
                pt2pt_mode=1;
            }
        }

        if (strcmp(argv[x],"-s")==0)
        {
            /* Read city file(s) */

            z=x+1;

            while (z<argc && argv[z][0] && argv[z][0]!='-' && cities<5)
            {
                strncpy(city_file[cities],argv[z],MAX_PATH_LEN-1);
                cities++;
                z++;
            }

            z--;
        }

        if (strcmp(argv[x],"-b")==0)
        {
            /* Read Boundary File(s) */

            z=x+1;

            while (z<argc && argv[z][0] && argv[z][0]!='-' && bfs<5)
            {
                strncpy(boundary_file[bfs],argv[z],MAX_PATH_LEN-1);
                bfs++;
                z++;
            }

            z--;
        }

        if (strcmp(argv[x],"-f")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&forced_freq);

                if (forced_freq<20.0)
                    forced_freq=0.0;

                if (forced_freq>20.0e3)
                    forced_freq=20.0e3;
            }
        }

        if (strcmp(argv[x],"-erp")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
            {
                sscanf(argv[z],"%lf",&forced_erp);

                if (forced_erp<0.0)
                    forced_erp=-1.0;
            }			 
        }

        if (strcmp(argv[x],"-ano")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
                strncpy(ano_filename,argv[z],253);
        }

        if (strcmp(argv[x],"-ani")==0)
        {
            z=x+1;

            if (z<argc && argv[z][0] && argv[z][0]!='-')
                strncpy(ani_filename,argv[z],253);
        }
    }

    sprintf(header,"\n\t\t--==[ Welcome To %s%s v%s ]==--\n\n",
            SPLAT_NAME, (appmode==APPMODE_HD?" HD":""), SPLAT_VERSION);

    /* Perform some error checking on the arguments
       and switches parsed from the command-line.
       If an error is encountered, print a message
       and exit gracefully. */

    if (txsites==0)
    {
        fprintf(stderr,"\n%c*** ERROR: No transmitter site(s) specified!\n\n",7);
        return -1;
    }

    errcount = 0;
    for (unsigned int x=0; x<txsites; x++)
    {
        if (tx_site[x].lat==91.0 && tx_site[x].lon==361.0)
        {
            fprintf(stderr,"\n*** ERROR: Transmitter site #%d not found!",x+1);
            errcount++;
        }
    }

    if (errcount)
    {
        fprintf(stderr,"%c\n\n",7);
        return -1;
    }

    /* XXX TODO: Read maxpagesides from commandline.
     * Check to see if we can support the requested gridsize.
     * If not, either issue warning or ratchet down.
     * For now, just set this to the maximum. */

    ippd = AppArraySize(appmode);
    ppd=(double)ippd;	/* pixels per degree (double)  */
    dpp=1.0/ppd;		/* degrees per pixel */
    mpi=ippd-1;		/* maximum pixel index per degree */


    if ((coverage+LRmap+ani_filename[0])==0 && rx_site.lat==91.0 && rx_site.lon==361.0)
    {
        if (max_range!=0.0 && txsites!=0)
        {
            /* Plot topographic map of radius "max_range" */

            map=0;
            topomap=1;
        }

        else
        {
            fprintf(stderr,"\n%c*** ERROR: No receiver site found or specified!\n\n",7);
            exit (-1);
        }
    }

    /* No major errors were detected.  Whew!  :-) */


    if (InitDEMs() < 0) {
        fprintf(stderr, "Couldn't allocate Digital Elevation Model array!\n");
        return -1;
    }

    /* Adjust input parameters if -metric option is used */

    if (metric)
    {
        altitudeLR/=METERS_PER_FOOT;	/* meters --> feet */
        max_range/=KM_PER_MILE;		/* kilometers --> miles */
        altitude/=METERS_PER_FOOT;	/* meters --> feet */
        clutter/=METERS_PER_FOOT;	/* meters --> feet */
    }

#ifndef _WIN32
    /* read the home_sdf_path */
    env=getenv("HOME");
    snprintf(buf,MAX_PATH_LEN-2,"%s/.splat_path",env);
    fd=fopen(buf,"r");
    if (fd!=NULL)
    {
        fgets(buf,MAX_PATH_LEN-2,fd);

        /* Remove <CR> and/or <LF> from buf */
        int x;
        for (x=0; buf[x]!=13 && buf[x]!=10 && buf[x]!=0 && x<(MAX_PATH_LEN-2); x++);
        buf[x]=0;

        strcpy(home_sdf_path,buf);
        fclose(fd);

        /* Ensure it has a trailing slash */
        size_t len = strlen(home_sdf_path);
        if (len > 0 && home_sdf_path[len-1] != PATHSEP) {
            home_sdf_path[len-1] = PATHSEP;
            home_sdf_path[len] = '\0';
        }

    }
#endif


    fprintf(stdout,"%s",header);
    fflush(stdout);

    if (ani_filename[0])
    {
        LoadLRP(tx_site[0],0); /* Get ERP status */
        LoadANO(ani_filename);

        for (unsigned int x=0; x<txsites && x<max_txsites; x++)
            PlaceMarker(tx_site[x]);

        if (rxsite)
            PlaceMarker(rx_site);

        if (bfs)
        {
            for (unsigned int x=0; x<bfs; x++)
                LoadBoundaries(boundary_file[x]);

            fprintf(stdout,"\n");
            fflush(stdout);
        }

        if (cities)
        {
            for (unsigned int x=0; x<cities; x++)
                LoadCities(city_file[x]);

            fprintf(stdout,"\n");
            fflush(stdout);
        }

        if (LR.erp==0.0)
            WriteImageLR(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);
        else
        {
            if (dbm)
                WriteImageDBM(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);
            else
                WriteImageSS(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);
        }

        FreeDEMs();
        return 0;
    }

    min_lat=90;
    max_lat=-90;

    min_lon=(int)floor(tx_site[0].lon);
    max_lon=(int)floor(tx_site[0].lon);

    for (unsigned int z=0; z<txsites && z<max_txsites; z++)
    {
        txlat=(int)floor(tx_site[z].lat);
        txlon=(int)floor(tx_site[z].lon);

        if (txlat<min_lat)
            min_lat=txlat;

        if (txlat>max_lat)
            max_lat=txlat;

        if (LonDiff(txlon,min_lon)<0.0)
            min_lon=txlon;

        if (LonDiff(txlon,max_lon)>=0.0)
            max_lon=txlon;
    }

    if (rxsite)
    {
        rxlat=(int)floor(rx_site.lat);
        rxlon=(int)floor(rx_site.lon);

        if (rxlat<min_lat)
            min_lat=rxlat;

        if (rxlat>max_lat)
            max_lat=rxlat;

        if (LonDiff(rxlon,min_lon)<0.0)
            min_lon=rxlon;

        if (LonDiff(rxlon,max_lon)>=0.0)
            max_lon=rxlon;
    }

    /* Load the required SDF files */ 

    LoadTopoData(max_lon, min_lon, max_lat, min_lat);

    if (area_mode || topomap)
    {
        for (unsigned int z=0; z<txsites && z<max_txsites; z++)
        {
            /* "Ball park" estimates used to load any additional
               SDF files required to conduct this analysis. */

            tx_range=sqrt(1.5*(tx_site[z].alt+GetElevation(tx_site[z])));

            if (LRmap)
                rx_range=sqrt(1.5*altitudeLR);
            else
                rx_range=sqrt(1.5*altitude);

            /* deg_range determines the maximum
               amount of topo data we read */

            deg_range=(tx_range+rx_range)/57.0;

            /* max_range regulates the size of the
               analysis.  A small, non-zero amount can
               be used to shrink the size of the analysis
               and limit the amount of topo data read by
               SPLAT!  A large number will increase the
               width of the analysis and the size of
               the map. */

            if (max_range==0.0)
                max_range=tx_range+rx_range;

            deg_range=max_range/57.0;

            /* Prevent the demand for a really wide coverage
               from allocating more "pages" than are available
               in memory. */
            deg_limit = .25;
            if (maxpagesides > 1) {
                deg_limit = 0.5 * (maxpagesides - 1.0);
            }

            if (fabs(tx_site[z].lat)<70.0)
                deg_range_lon=deg_range/cos(DEG2RAD*tx_site[z].lat);
            else
                deg_range_lon=deg_range/cos(DEG2RAD*70.0);

            /* Correct for squares in degrees not being square in miles */  

            if (deg_range>deg_limit)
                deg_range=deg_limit;

            if (deg_range_lon>deg_limit)
                deg_range_lon=deg_limit;

            north_min=(int)floor(tx_site[z].lat-deg_range);
            north_max=(int)floor(tx_site[z].lat+deg_range);

            west_min=(int)floor(tx_site[z].lon-deg_range_lon);

            while (west_min<0)
                west_min+=360;

            while (west_min>=360)
                west_min-=360;

            west_max=(int)floor(tx_site[z].lon+deg_range_lon);

            while (west_max<0)
                west_max+=360;

            while (west_max>=360)
                west_max-=360;

            if (north_min<min_lat)
                min_lat=north_min;

            if (north_max>max_lat)
                max_lat=north_max;

            if (LonDiff(west_min,min_lon)<0.0)
                min_lon=west_min;

            if (LonDiff(west_max,max_lon)>=0.0)
                max_lon=west_max;
        }

        /* Load any additional SDF files, if required */ 

        LoadTopoData(max_lon, min_lon, max_lat, min_lat);
    }

#ifndef _WIN32
    if (udt_file[0])
        LoadUDT(udt_file);
#endif

    /**** Set terrain elevation to zero for sites providing AMSL antenna heights ****

      for (z=0; z<txsites && z<max_txsites; z++)
      {
      if (tx_site[z].amsl_flag)
      x=SetElevation(tx_site[z].lat,tx_site[z].lon,0.0);
      }

      if (rxsite && rx_site.amsl_flag)
      SetElevation(rx_site.lat,rx_site.lon,0.0);

     *********************************************************************************/

    /***** Let the SPLATting begin! *****/

    if (pt2pt_mode)
    {
        char longley_ext[4], terrain_ext[4], elevation_ext[4], height_ext[4], norm_height_ext[4];

        PlaceMarker(rx_site);

        if (longley_plot) {
            stripExtension(longley_file, longley_ext, 4);
            if (strlen(longley_ext) == 0) {
                strcpy(longley_ext, "png");
            }
        }

        if (terrain_plot) {
            stripExtension(terrain_file, terrain_ext, 4);
            if (strlen(terrain_ext) == 0) {
                strcpy(terrain_ext, "png");
            }
        }

        if (elevation_plot) {
            stripExtension(elevation_file, elevation_ext, 4);
            if (strlen(elevation_ext) == 0) {
                strcpy(elevation_ext, "png");
            }
        }

        if (height_plot) {
            stripExtension(height_file, height_ext, 4);
            if (strlen(height_ext) == 0) {
                strcpy(height_ext, "png");
            }
        }

        if (norm_height_plot) {
            stripExtension(norm_height_file, norm_height_ext, 4);
            if (strlen(norm_height_ext) == 0) {
                strcpy(norm_height_ext, "png");
            }
        }

        for (unsigned int x=0; x<txsites && x<4; x++)
        {
            PlaceMarker(tx_site[x]);

            if (nolospath==0)
            {
                switch (x)
                {
                    case 0:
                        PlotPath(tx_site[x],rx_site,1);
                        break;

                    case 1:
                        PlotPath(tx_site[x],rx_site,8);
                        break;

                    case 2:
                        PlotPath(tx_site[x],rx_site,16);
                        break;

                    case 3:
                        PlotPath(tx_site[x],rx_site,32);
                }
            }

            if (nositereports==0)
                SiteReport(tx_site[x]);

            if (kml)
                WriteKML(tx_site[x],rx_site);

            if (nositereports==0)
            {
                if (txsites>1)
                    snprintf(buf,MAX_PATH_LEN,"%s-%d.%s", longley_file, x, longley_ext);
                else
                    snprintf(buf,MAX_PATH_LEN,"%s.%s", longley_file, longley_ext);

                LoadLRP(tx_site[x],longley_plot);
                PathReport(tx_site[x],rx_site,buf,longley_plot);
            }

            if (terrain_plot)
            {
                if (txsites>1)
                    snprintf(buf,MAX_PATH_LEN,"%s-%d.%s",terrain_file,x, terrain_ext);
                else
                    snprintf(buf,MAX_PATH_LEN,"%s.%s",terrain_file,terrain_ext);

                GraphTerrain(tx_site[x],rx_site,buf);
            }

            if (elevation_plot)
            {
                if (txsites>1)
                    snprintf(buf,MAX_PATH_LEN,"%s-%d.%s",elevation_file,x,elevation_ext);
                else
                    snprintf(buf,MAX_PATH_LEN,"%s.%s",elevation_file,elevation_ext);

                GraphElevation(tx_site[x],rx_site,buf);
            }

            if (height_plot)
            {
                if (txsites>1)
                    snprintf(buf,MAX_PATH_LEN,"%s-%d.%s",height_file,x,height_ext);
                else
                    snprintf(buf,MAX_PATH_LEN,"%s.%s",height_file,height_ext);

                GraphHeight(tx_site[x],rx_site,buf,fresnel_plot,0);
            }

            if (norm_height_plot)
            {
                if (txsites>1)
                    snprintf(buf,MAX_PATH_LEN,"%s-%d.%s",norm_height_file,x,norm_height_ext);
                else
                    snprintf(buf,MAX_PATH_LEN,"%s.%s",norm_height_file,norm_height_ext);

                GraphHeight(tx_site[x],rx_site,buf,fresnel_plot,1);
            }
        }
    }

    if (area_mode && topomap==0)
    {
        for (unsigned int x=0; x<txsites && x<max_txsites; x++)
        {
            if (coverage)
                PlotLOSMap(tx_site[x],altitude, multithread);

            else if (LoadLRP(tx_site[x],1))
                PlotLRMap(tx_site[x],altitudeLR,ano_filename, multithread);

            SiteReport(tx_site[x]);
        }
    }

    if (map || topomap)
    {
        /* Label the map */

        if (kml==0)
        {
            for (unsigned int x=0; x<txsites && x<max_txsites; x++)
                PlaceMarker(tx_site[x]);
        }

        if (cities)
        {

            for (unsigned int y=0; y<cities; y++)
                LoadCities(city_file[y]);

            fprintf(stdout,"\n");
            fflush(stdout);
        }

        /* Load city and county boundary data files */

        if (bfs)
        {
            for (unsigned int y=0; y<bfs; y++)
                LoadBoundaries(boundary_file[y]);

            fprintf(stdout,"\n");
            fflush(stdout);
        }

        /* Plot the map */

        if (coverage || pt2pt_mode || topomap)
            WriteImage(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);

        else
        {
            if (LR.erp==0.0)
                WriteImageLR(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);
            else
                if (dbm)
                    WriteImageDBM(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);
                else
                    WriteImageSS(mapfile,imagetype,geo,kml,ngs,tx_site,txsites);
        }
    }

    if (command_line_log && strlen(logfile)>0)
    {
        fd=fopen(logfile,"w");

        if (fd!=NULL)
        {
            for (int x=0; x<argc; x++)
                fprintf(fd,"%s ",argv[x]);

            fprintf(fd,"\n");

            fclose(fd);

            fprintf(stdout,"\nCommand-line parameter log written to: \"%s\"\n",logfile);

        }
    }

    printf("\n");

    FreeDEMs();

    /* That's all, folks! */

    return 0;
}
