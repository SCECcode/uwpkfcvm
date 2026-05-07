/**
 * @file uwpkfcvm.h
 * @brief Main header file for uwpkfcvm library.
 * @version 1.0
 *
 * Delivers the uwpkfcvm model 
 * base on original uwpkfcvm
 *
 */
#ifndef UWPKFCVM_H
#define UWPKFCVM_H

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "kdtree_util.h"
#include "proj.h"

// Constants
#ifndef M_PI
	/** Defines pi */
	#define M_PI 3.14159265358979323846
#endif
#define DEG_TO_RAD M_PI / 180.0

/** Defines a return value of success */
#define SUCCESS 0
/** Defines a return value of failure */
#define FAIL 1

/* config string */
#define UWPKFCVM_CONFIG_MAX 1000

// Structures
/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct uwpkfcvm_point_t {
	/** Longitude member of the point */
	double longitude;
	/** Latitude member of the point */
	double latitude;
	/** Depth member of the point */
	double depth;
} uwpkfcvm_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct uwpkfcvm_properties_t {
	/** P-wave velocity in meters per second */
	double vp;
	/** S-wave velocity in meters per second */
	double vs;
	/** Density in g/m^3 */
	double rho;
	/** Qp */
	double qp;
	/** Qs */
	double qs;
} uwpkfcvm_properties_t;

/** The CVM-S5 configuration structure. */
typedef struct uwpkfcvm_configuration_t {
	/** The zone of UTM projection */
	int utm_zone;
	/** The model directory */
	char model_dir[128];
	/** Number of x points */
	int nx;
	/** Number of y points */
	int ny;
	/** Number of z points */
	int nz;
	/** Depth in meters */
	double depth;
	/** Top left corner easting in UTM projection */
	double top_left_corner_e;
	/** Top left corner northing in UTM projection */
	double top_left_corner_n;
	/** Top right corner easting in UTM projection */
	double top_right_corner_e;
	/** Top right corner northing in UTM projection */
	double top_right_corner_n;
	/** Bottom left corner easting in UTM projection */
	double bottom_left_corner_e;
	/** Bottom left corner northing in UTM projection */
	double bottom_left_corner_n;
	/** Bottom right corner easting in UTM projection */
	double bottom_right_corner_e;
	/** Bottom right corner northing in UTM projection */
	double bottom_right_corner_n;
        /** The data access seek method, fast-X, or fast-Y */
        char seek_axis[128];
        /** The data seek direction, bottom-up, or top-down */
        char seek_direction[128];
        /** trilinear interploation; */
        int interpolation;
        /** add 1d model using nearest boundary points */
        int 1d_background;
} uwpkfcvm_configuration_t;

/** The model structure which points to available portions of the model. */
typedef struct uwpkfcvm_model_t {
	KDlld *pnts;       // raw model data
	int *pnts_zero_depth; // track which pnts has depth/z = 0
			      
        KDVec3 *v3pnts;    // x,y,z
        KDNode3 *v3nodes;  // tree v3nodes for vec3

        KDVec2 *v2pnts; // utm_e, utm_n, just 1 layer 
        KDVec2 *v2pnts_boundary;
        int boundary_size;

        KDVec2 *v2hull;

	int v2hull_size;
        int zero_depth_cnt; // should be NX * NY			


} uwpkfcvm_model_t;

// UCVM API Required Functions
#ifdef DYNAMIC_LIBRARY

/** Initializes the model */
int model_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int model_finalize();
/** Returns version information */
int model_version(char *ver, int len);
/** Queries the model */
int model_query(uwpkfcvm_point_t *points, uwpkfcvm_properties_t *data, int numpts);
int model_config(char **config, int *sz);

int (*get_model_init())(const char *, const char *);
int (*get_model_query())(uwpkfcvm_point_t *, uwpkfcvm_properties_t *, int);
int (*get_model_finalize())();
int (*get_model_version())(char *, int);
int (*get_model_config())(char **, int*);

#endif

// uwpkfcvm Related Functions

/** Initializes the model */
int uwpkfcvm_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int uwpkfcvm_finalize();
/** Returns version information */
int uwpkfcvm_version(char *ver, int len);
/** Queries the model */
int uwpkfcvm_query(uwpkfcvm_point_t *points, uwpkfcvm_properties_t *data, int numpts);

void uwpkfcvm_setdebug();

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int uwpkfcvm_read_configuration(char *file, uwpkfcvm_configuration_t *config);
int uwpkfcvm_dump_configuration(uwpkfcvm_configuration_t *config);
/** Prints out the error string. */
void uwpkfcvm_print_error(char *err);
/** Retrieves the value at a specified grid point in the model. */
void uwpkfcvm_read_properties2(int x, int y, int z, uwpkfcvm_properties_t *data);
void uwpkfcvm_read_properties(uwpkfcvm_model_t *model, int index, uwpkfcvm_properties_t *data);
/** Attempts to malloc the model size in memory and read it in. */
int uwpkfcvm_reading_model(uwpkfcvm_model_t *model);

// Interpolation Functions
/** Linearly interpolates two uwpkfcvm_properties_t structures */
void uwpkfcvm_linear_interpolation(double percent, uwpkfcvm_properties_t *x0, uwpkfcvm_properties_t *x1, uwpkfcvm_properties_t *ret_properties);
/** Bilinearly interpolates the properties. */
void uwpkfcvm_bilinear_interpolation(double x_percent, double y_percent, uwpkfcvm_properties_t *four_points, uwpkfcvm_properties_t *ret_properties);
/** Trilinearly interpolates the properties. */
void uwpkfcvm_trilinear_interpolation(double x_percent, double y_percent, double z_percent, uwpkfcvm_properties_t *eight_points,
							 uwpkfcvm_properties_t *ret_properties);

// from uwpkfcvm_util.c
void setup_model(uwpkfcvm_model_t *model, int cnt);
void free_model(uwpkfcvm_model_t *model);
void load_model(uwpkfcvm_model_t *model, int NX, int NY, int NZ, FILE *fp);
int in_model(uwpkfcvm_model_t *model, double lat, double lon, double depth);
int nearest_neighbor(uwpkfcvm_model_t *model, double lat, double lon, double depth);
double vs_by_location(uwpkfcvm_model_t *model,int loc);
double vp_by_location(uwpkfcvm_model_t *model,int loc);

#endif // UWPKFCVM_H
