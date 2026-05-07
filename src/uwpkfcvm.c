/**
 * @file uwpkfcvm.c
 * @brief Main file for uwpkfcvm model
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers the uwpkfcvm model which base from 
 * original linthurber dataset
 * origin: lower left, fast-x, top-down 
 */

#include "limits.h"
#include "ucvm_model_dtypes.h"
#include "uwpkfcvm.h"

#include <assert.h>

int uwpkfcvm_ucvm_debug_detail=0;
int uwpkfcvm_ucvm_debug=0;
FILE *stderrfp=NULL;

/** The config of the model */
char *uwpkfcvm_config_string=NULL;
int uwpkfcvm_config_sz=0;

// Constants
/** The version of the model. */
const char *uwpkfcvm_version_string = "uwpkfcvm";

// Variables
/** Set to 1 when the model is ready for query. */
int uwpkfcvm_is_initialized = 0;

char uwpkfcvm_data_directory[128];

/** Configuration parameters. */
uwpkfcvm_configuration_t *uwpkfcvm_configuration;
/** Holds pointers to the velocity model data and its derived data */
uwpkfcvm_model_t *uwpkfcvm_velocity_model;

/** Proj coordinate transformation objects. can go from geo <-> utm */
PJ *uwpkfcvm_geo2utm = NULL;

/** The cosine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double uwpkfcvm_cos_rotation_angle = 0;
/** The sine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double uwpkfcvm_sin_rotation_angle = 0;

/** The height of this model's region, in meters. */
double uwpkfcvm_total_height_m = 0;
/** The width of this model's region, in meters. */
double uwpkfcvm_total_width_m = 0;
/**
 * Initializes the uwpkfcvm plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int uwpkfcvm_init(const char *dir, const char *label) {
    int tempVal = 0;
    char configbuf[512];
    double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

    if(uwpkfcvm_ucvm_debug) {
      stderrfp = fopen("uwpkfcvm_debug.log", "w+");
      fprintf(stderrfp," ===== START ===== \n");
    }


    // Initialize variables.
    uwpkfcvm_configuration = calloc(1, sizeof(uwpkfcvm_configuration_t));
    uwpkfcvm_velocity_model = calloc(1, sizeof(uwpkfcvm_model_t));

    uwpkfcvm_config_string = calloc(UWPKFCVM_CONFIG_MAX, sizeof(char));
    uwpkfcvm_config_string[0]='\0';
    uwpkfcvm_config_sz=0;

    // Configuration file location.
    sprintf(configbuf, "%s/model/%s/data/config", dir, label);

    // Read the uwpkfcvm_configuration file.
    if (uwpkfcvm_read_configuration(configbuf, uwpkfcvm_configuration) != SUCCESS)
        return FAIL;

    // Set up the data directory.
    sprintf(uwpkfcvm_data_directory, "%s/model/%s/data/%s/", dir, label, uwpkfcvm_configuration->model_dir);

    // We need to convert the point from lat, lon to UTM, let's set it up.
    char uwpkfcvm_projstr[64];
    snprintf(uwpkfcvm_projstr, 64, "+proj=utm +ellps=clrk66 +zone=%d +datum=NAD27 +units=m +no_defs", uwpkfcvm_configuration->utm_zone);
    if (!(uwpkfcvm_geo2utm = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", uwpkfcvm_projstr, NULL))) {
        uwpkfcvm_print_error("Could not set up Proj transformation from EPSG:4326 to UTM.");
        uwpkfcvm_print_error((char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (UCVM_CODE_ERROR);
    }

    // Can we allocate the model, or parts of it, to memory. If so, we do.
    tempVal = uwpkfcvm_reading_model(uwpkfcvm_velocity_model);

    if (tempVal == SUCCESS) {
//        fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
//        fprintf(stderr, "hard disk may result in slow performance.\n");
    } else if (tempVal == FAIL) {
        uwpkfcvm_print_error("No model file was found to read from.");
        return FAIL;
    }

    // setup config_string 
    sprintf(uwpkfcvm_config_string,"config = %s\n",configbuf);
    uwpkfcvm_config_sz=1;


    // Let everyone know that we are initialized and ready for business.
    uwpkfcvm_is_initialized = 1;

    return SUCCESS;
}

/**
 * Queries uwpkfcvm at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, rho, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int uwpkfcvm_query(uwpkfcvm_point_t *points, uwpkfcvm_properties_t *data, int numpoints) {
    int i = 0;

    double point_u = 0, point_v = 0;
    double point_x = 0, point_y = 0; 
				   
    int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
    double x_percent = 0, y_percent = 0, z_percent = 0;

    uwpkfcvm_properties_t surrounding_points[8];
    int zone = uwpkfcvm_configuration->utm_zone;

    for (i = 0; i < numpoints; i++) {

        // We need to be below the surface to service this query.
        if (points[i].depth < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

	// is it in model ??
        if( ! in_model(uwpkfcvm_velocity_model, points[i].latitude, points[i].longitude, points[i].depth) ) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
            } else { 
	    // query for nearest data point
                int modelindex=nearest_neighbor(uwpkfcvm_velocity_model, points[i].latitude, points[i].longitude, points[i].depth);
                //if(uwpkfcvm_configuration->interpolation) { // do nothing for now }
                uwpkfcvm_read_properties(uwpkfcvm_velocity_model, modelindex, &(data[i]));    
        }
    }

    return SUCCESS;
}


/**
 * Calculates the vs based off of Vp. Base on Brocher's formulae
 *
 * https://pubs.usgs.gov/of/2005/1317/of2005-1317.pdf
 *
 * @param vp
 * @return Vs, in km.
 * Vs derived from Vp, Brocher (2005) eqn 1.
 * [eqn. 1] Vs (km/s) = 0.7858 – 1.2344Vp + 0.7949Vp2 – 0.1238Vp3 + 0.0064Vp4.
 * Equation 1 is valid for 1.5 < Vp < 8 km/s.
 */
double uwpkfcvm_calculate_vs(double vp) {
     double retVal ;

     vp = vp * 0.001;
     double t1= (vp * 1.2344);
     double t2= ((vp * vp)* 0.7949);
     double t3= ((vp * vp * vp) * 0.1238);
     double t4= ((vp * vp * vp * vp) * 0.0064);
     retVal = 0.7858 - t1 + t2 - t3 + t4;
     retVal = retVal * 1000.0;

     return retVal;
}

/**
 * Calculates the density based off of Vp. Base on Brocher's formulae
 *
 * @param vp
 * @return Density, in g/m^3.
 * [eqn. 6] r (g/cm3) = 1.6612Vp – 0.4721Vp2 + 0.0671Vp3 – 0.0043Vp4 + 0.000106Vp5.
 * Equation 6 is the “Nafe-Drake curve” (Ludwig et al., 1970).
 * start with vp in km
 */
double uwpkfcvm_calculate_density(double vp) {
     double retVal ;

     vp = vp * 0.001;
     double t1 = (vp * 1.6612);
     double t2 = ((vp * vp ) * 0.4721);
     double t3 = ((vp * vp * vp) * 0.0671);
     double t4 = ((vp * vp * vp * vp) * 0.0043);
     double t5 = ((vp * vp * vp * vp * vp) * 0.000106);
     retVal = t1 - t2 + t3 - t4 + t5;
     if (retVal < 1.0) {
       retVal = 1.0;
     }
     retVal = retVal * 1000.0;
     return retVal;
}

void uwpkfcvm_read_properties(uwpkfcvm_model_t *model, int location, uwpkfcvm_properties_t *data) {  

    double vs=vs_by_location(model,location);

    data->vs = vs;
    data->vp=vp_by_location(model,location);
    /* Calculate density */
    if (data->vp > 0.0) { data->rho=uwpkfcvm_calculate_density(data->vp); }
}

/**
 */
void uwpkfcvm_setdebug() {
   uwpkfcvm_ucvm_debug=1;
}


/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int uwpkfcvm_finalize() {

    proj_destroy(uwpkfcvm_geo2utm);
    uwpkfcvm_geo2utm = NULL;

    if (uwpkfcvm_velocity_model) {
      free_model(uwpkfcvm_velocity_model);
      free(uwpkfcvm_velocity_model);
    }
    if (uwpkfcvm_configuration) free(uwpkfcvm_configuration);

    return SUCCESS;

    if(uwpkfcvm_ucvm_debug) {
      fprintf(stderrfp,"DONE:\n");
      fclose(stderrfp);
    }

}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int uwpkfcvm_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(uwpkfcvm_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, uwpkfcvm_version_string, verlen);
  return 0;
}

/**
 * Returns the model config information.
 *
 * @param key Config key string to return.
 * @param sz number of config terms.
 * @return Zero
 */
int uwpkfcvm_config(char **config, int *sz)
{
  int len=strlen(uwpkfcvm_config_string);
  if(len > 0) {
    *config=uwpkfcvm_config_string;
    *sz=uwpkfcvm_config_sz;
    return SUCCESS;
  }
  return FAIL;
}


/**
 * Reads the uwpkfcvm_configuration file describing the various properties of CVM-S5 and populates
 * the uwpkfcvm_configuration struct. This assumes uwpkfcvm_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The uwpkfcvm_configuration file location on disk to read.
 * @param config The uwpkfcvm_configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int uwpkfcvm_read_configuration(char *file, uwpkfcvm_configuration_t *config) {
    FILE *fp = fopen(file, "r");
    char key[40];
    char value[80];
    char line_holder[128];

    // If our file pointer is null, an error has occurred. Return fail.
    if (fp == NULL) {
        uwpkfcvm_print_error("Could not open the uwpkfcvm_configuration file.");
        return FAIL;
    }

    // Read the lines in the uwpkfcvm_configuration file.
    while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
        if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
            sscanf(line_holder, "%s = %s", key, value);

            // Which variable are we editing?
            if (strcmp(key, "utm_zone") == 0) config->utm_zone = atoi(value);
            if (strcmp(key, "model_dir") == 0) sprintf(config->model_dir, "%s", value);
            if (strcmp(key, "nx") == 0) config->nx = atoi(value);
            if (strcmp(key, "ny") == 0) config->ny = atoi(value);
            if (strcmp(key, "nz") == 0) config->nz = atoi(value);
            if (strcmp(key, "depth") == 0) config->depth = atof(value);
            if (strcmp(key, "top_left_corner_e") == 0) config->top_left_corner_e = atof(value);
            if (strcmp(key, "top_left_corner_n") == 0) config->top_left_corner_n = atof(value);
            if (strcmp(key, "top_right_corner_e") == 0) config->top_right_corner_e = atof(value);
            if (strcmp(key, "top_right_corner_n") == 0) config->top_right_corner_n = atof(value);
            if (strcmp(key, "bottom_left_corner_e") == 0) config->bottom_left_corner_e = atof(value);
            if (strcmp(key, "bottom_left_corner_n") == 0) config->bottom_left_corner_n = atof(value);
            if (strcmp(key, "bottom_right_corner_e") == 0) config->bottom_right_corner_e = atof(value);
            if (strcmp(key, "bottom_right_corner_n") == 0) config->bottom_right_corner_n = atof(value);
            if (strcmp(key, "seek_axis") == 0) sprintf(config->seek_axis, "%s", value);
            if (strcmp(key, "seek_direction") == 0) sprintf(config->seek_direction, "%s", value);
            if (strcmp(key, "interpolation") == 0) { 
                config->interpolation=0;
                if (strcmp(value,"on") == 0) config->interpolation=1;
            }
            if (strcmp(key, "1d_background") == 0) { 
                config->id_background=0;
                if (strcmp(value,"on") == 0) config->id_background=1;
            }
        }
    }

    // Have we set up all uwpkfcvm_configuration parameters?
    if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' || 
        config->seek_direction[0] == '\0' || config->seek_axis[0] == '\0' ||
        config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
        config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
        config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ) {
        uwpkfcvm_print_error("One of uwpkfcvm_configuration parameter not specified. Please check your uwpkfcvm_configuration file.");
        return FAIL;
    }

    fclose(fp);

    return SUCCESS;
}

/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void uwpkfcvm_print_error(char *err) {
    fprintf(stderr, "An error has occurred while executing uwpkfcvm. The error was: %s\n",err);
    fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
    fprintf(stderr, "about the computer you are running uwpkfcvm on (Linux, Mac, etc.).\n");
}

/**
 * Check if the data is too big to be loaded internally (exceed maximum
 * allowable by a INT variable)
 *
 */
static int too_big() {
        long max_size= (long) (uwpkfcvm_configuration->nx) * uwpkfcvm_configuration->ny * uwpkfcvm_configuration->nz;
        long delta= max_size - INT_MAX;

    if( delta > 0) {
        return 1;
        } else {
        return 0;
        }
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return FAIL 1, SUCCESS if processed okay, 0
 * is not in memory, FAIL if no file found.
 */
int uwpkfcvm_reading_model(uwpkfcvm_model_t *model) {
	
    int file_count = 0;
    char current_file[128];
    FILE *fp;

    // Let's see what data we actually have.
    sprintf(current_file, "%s/parkfield.txt", uwpkfcvm_data_directory);
    if (access(current_file, R_OK) == 0) {
       fp = fopen(current_file, "rb");
       load_model(model, uwpkfcvm_configuration->nx,
		       uwpkfcvm_configuration->ny, uwpkfcvm_configuration->nz, fp);
       return SUCCESS;
    }
    return FAIL;
}

// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls uwpkfcvm_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
    return uwpkfcvm_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls uwpkfcvm_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(uwpkfcvm_point_t *points, uwpkfcvm_properties_t *data, int numpoints) {
    return uwpkfcvm_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls uwpkfcvm_finalize.
 *
 * @return Success
 */
int model_finalize() {
    return uwpkfcvm_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls uwpkfcvm_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
    return uwpkfcvm_version(ver, len);
}

/**
 * Version function loaded and called by the UCVM library. Calls uwpkfcvm_config.
 *
 * @param config Config string to return.
 * @param sz number of config terms
 * @return Zero
 */
int model_config(char **config, int *sz) {
    return uwpkfcvm_config(config, sz);
}


int (*get_model_init())(const char *, const char *) {
        return &uwpkfcvm_init;
}
int (*get_model_query())(uwpkfcvm_point_t *, uwpkfcvm_properties_t *, int) {
         return &uwpkfcvm_query;
}
int (*get_model_finalize())() {
         return &uwpkfcvm_finalize;
}
int (*get_model_version())(char *, int) {
         return &uwpkfcvm_version;
}
int (*get_model_config())(char **, int*) {
    return &uwpkfcvm_config;
}



#endif
