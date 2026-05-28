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
/** Holds pointers to the velocity model data OR indicates it can be read from file. */
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

    // Can we allocate the model, or parts of it, to memory. If so, we do.
    tempVal = uwpkfcvm_try_reading_model(uwpkfcvm_velocity_model);

    if (tempVal == SUCCESS) {
//        fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
//        fprintf(stderr, "hard disk may result in slow performance.\n");
    } else if (tempVal == FAIL) {
        uwpkfcvm_print_error("No model file was found to read from.");
        return FAIL;
    }

    // We need to convert the point from lat, lon to UTM, let's set it up.
    char uwpkfcvm_projstr[64];
    snprintf(uwpkfcvm_projstr, 64, "+proj=utm +ellps=clrk66 +zone=%d +datum=NAD27 +units=m +no_defs", uwpkfcvm_configuration->utm_zone);
    if (!(uwpkfcvm_geo2utm = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", uwpkfcvm_projstr, NULL))) {
        uwpkfcvm_print_error("Could not set up Proj transformation from EPSG:4326 to UTM.");
        uwpkfcvm_print_error((char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (UCVM_CODE_ERROR);
    }


    // In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
    // corner is at (0m,0m). Our box's height is total_height_m and total_width_m. We then rotate the
    // point so that is is somewhere between (0,0) and (total_width_m, total_height_m). How far along
    // the X and Y axis determines which grid points we use for the interpolation routine.

    // Calculate the rotation angle of the box.
    north_height_m = uwpkfcvm_configuration->top_left_corner_n - uwpkfcvm_configuration->bottom_left_corner_n;
    east_width_m = uwpkfcvm_configuration->top_left_corner_e - uwpkfcvm_configuration->bottom_left_corner_e;

    // Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
    rotation_angle = atan(east_width_m / north_height_m);

    uwpkfcvm_cos_rotation_angle = cos(rotation_angle);
    uwpkfcvm_sin_rotation_angle = sin(rotation_angle);

    uwpkfcvm_total_height_m = sqrt(pow(uwpkfcvm_configuration->top_left_corner_n - uwpkfcvm_configuration->bottom_left_corner_n, 2.0f) +
          pow(uwpkfcvm_configuration->top_left_corner_e - uwpkfcvm_configuration->bottom_left_corner_e, 2.0f));
    uwpkfcvm_total_width_m  = sqrt(pow(uwpkfcvm_configuration->top_right_corner_n - uwpkfcvm_configuration->top_left_corner_n, 2.0f) +
          pow(uwpkfcvm_configuration->top_right_corner_e - uwpkfcvm_configuration->top_left_corner_e, 2.0f));

    if(uwpkfcvm_ucvm_debug) {
      fprintf(stderrfp,"north_height %lf east_width %lf\n", north_height_m, east_width_m);
      fprintf(stderrfp,"totol height %lf total width %lf\n", uwpkfcvm_total_height_m, uwpkfcvm_total_width_m);
      fprintf(stderrfp,"cos angle %lf sin angle %lf\n", uwpkfcvm_cos_rotation_angle, uwpkfcvm_sin_rotation_angle);
    }

    // setup config_string 
    sprintf(uwpkfcvm_config_string,"config = %s\n",configbuf);
    uwpkfcvm_config_sz=1;


    // Let everyone know that we are initialized and ready for business.
    uwpkfcvm_is_initialized = 1;

    return SUCCESS;
}

static int to_utm(double lon, double lat, double *point_u, double *point_v) {
    PJ_COORD xyzSrc = proj_coord(lat, lon, 0.0, HUGE_VAL);
    PJ_COORD xyzDest = proj_trans(uwpkfcvm_geo2utm, PJ_FWD, xyzSrc);
    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming latitude=%.4f, longitude=%.4f to UTM.\n",
              lat, lon);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return UCVM_CODE_ERROR;
    }
    *point_u = xyzDest.xyzt.x;
    *point_v = xyzDest.xyzt.y;
    return err;
}

static int to_geo(double point_u, double point_v, double *lon, double *lat) {
    PJ_COORD xyzSrc;
    xyzSrc.xyzt.x=point_u;
    xyzSrc.xyzt.y=point_v;
    PJ_COORD xyzDest = proj_trans(uwpkfcvm_geo2utm, PJ_INV, xyzSrc);
    
    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming u=%.4f, v=%.4f to Geo.\n",
              point_u, point_v);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return UCVM_CODE_ERROR;
    }
    *lon=xyzDest.lp.lam;
    *lat=xyzDest.lp.phi;
    return err;
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

	// lon,lat,u,v			     
	to_utm(points[i].longitude, points[i].latitude, &point_u, &point_v);

if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   left_e %lf left_n %lf\n", 
                      uwpkfcvm_configuration->bottom_left_corner_e, uwpkfcvm_configuration->bottom_left_corner_n); }

if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   lon %lf lat %lf\n", points[i].longitude, points[i].latitude); }
if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   point_u %lf point_v %lf\n", point_u, point_v); }

        // Point within rectangle.
        point_u -= uwpkfcvm_configuration->bottom_left_corner_e;
        point_v -= uwpkfcvm_configuration->bottom_left_corner_n;
if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"2  point_u %lf point_v %lf\n", point_u, point_v); }

        // We need to rotate that point, the number of degrees we calculated above.
        point_x = uwpkfcvm_cos_rotation_angle * point_u - uwpkfcvm_sin_rotation_angle * point_v;
        point_y = uwpkfcvm_sin_rotation_angle * point_u + uwpkfcvm_cos_rotation_angle * point_v;
if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   point_x %lf point_y  %lf\n", point_x, point_y); }

        // Which point base point does that correspond to?
        load_x_coord = round(point_x / uwpkfcvm_total_width_m * (uwpkfcvm_configuration->nx - 1));

        load_y_coord = round(point_y / uwpkfcvm_total_height_m * (uwpkfcvm_configuration->ny - 1));

        // And on the Z-axis? XXX
	
        load_z_coord = (uwpkfcvm_configuration->depth / uwpkfcvm_configuration->depth_interval) -
                       round(points[i].depth / uwpkfcvm_configuration->depth_interval);

if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   load_x_coord %d load_y_coord %d load_z_coord %d\n", load_x_coord,load_y_coord,load_z_coord); }

        // Are we outside the model's X and Y boundaries?
	// and also outside of z
        if (load_x_coord > uwpkfcvm_configuration->nx - 1 || load_y_coord > uwpkfcvm_configuration->ny - 1 || load_x_coord < 0 || load_y_coord < 0 | load_z_coord < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

        if(uwpkfcvm_configuration->interpolation) {

          // Get the X, Y, and Z percentages for the bilinear or trilinear interpolation below.
          double x_interval=(uwpkfcvm_configuration->nx > 1) ?
                     uwpkfcvm_total_width_m / (uwpkfcvm_configuration->nx-1):uwpkfcvm_total_width_m;
          double y_interval=(uwpkfcvm_configuration->ny > 1) ?
                     uwpkfcvm_total_height_m / (uwpkfcvm_configuration->ny-1):uwpkfcvm_total_height_m;

          x_percent = fmod(point_u, x_interval) / x_interval;
          y_percent = fmod(point_v, y_interval) / y_interval;
          z_percent = fmod(points[i].depth, uwpkfcvm_configuration->depth_interval) / uwpkfcvm_configuration->depth_interval;

          // Read all the surrounding point properties.
          uwpkfcvm_read_properties(load_x_coord, load_y_coord, load_z_coord, &(surrounding_points[0]));    // Orgin.
          uwpkfcvm_read_properties(load_x_coord + 1, load_y_coord, load_z_coord, &(surrounding_points[1]));    // Orgin + 1x
          uwpkfcvm_read_properties(load_x_coord, load_y_coord + 1, load_z_coord, &(surrounding_points[2]));    // Orgin + 1y
          uwpkfcvm_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord, &(surrounding_points[3]));    // Orgin + x + y, forms top plane.
          uwpkfcvm_read_properties(load_x_coord, load_y_coord, load_z_coord - 1, &(surrounding_points[4]));    // Bottom plane origin
          uwpkfcvm_read_properties(load_x_coord + 1, load_y_coord, load_z_coord - 1, &(surrounding_points[5]));    // +1x
          uwpkfcvm_read_properties(load_x_coord, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));    // +1y
          uwpkfcvm_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));    // +x +y, forms bottom plane.
  
          uwpkfcvm_trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
          } else {
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp,"direct call, no interpolation\n"); }
              uwpkfcvm_read_properties(load_x_coord, load_y_coord, load_z_coord, &(data[i]));    // Orgin.
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

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void uwpkfcvm_read_properties(int x, int y, int z, uwpkfcvm_properties_t *data) {
    // Set everything to -1 to indicate not found.
    data->vp = -1;
    data->vs = -1;
    data->rho = -1;
    data->qp = -1;
    data->qs = -1;

if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp,"read_properties index: x(%d) y(%d) z(%d)",x,y,z); }
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp," from (%d:%d:%d)\n",
	          uwpkfcvm_configuration->nx,uwpkfcvm_configuration->ny,uwpkfcvm_configuration->nz); }
    float *ptr = NULL;
    FILE *fp = NULL;
    long location = 0;


    // the z is inverted at line #145
    if ( strcmp(uwpkfcvm_configuration->seek_axis, "fast-y") == 0 ||
                 strcmp(uwpkfcvm_configuration->seek_axis, "fast-Y") == 0 ) { // fast-y,  uwpkfcvm 
        if(strcmp(uwpkfcvm_configuration->seek_direction, "bottom-up") == 0) { 
            location = ((long) z * uwpkfcvm_configuration->nx * uwpkfcvm_configuration->ny) + (x * uwpkfcvm_configuration->ny) + y;
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp," >>> LOCATION==%ld(fast-y, bottom-up)\n", location); }
            } else { // nz starts from 0 up to nz-1
                location = ((long)((uwpkfcvm_configuration->nz -1) - z) * uwpkfcvm_configuration->nx * uwpkfcvm_configuration->ny) + (x * uwpkfcvm_configuration->ny) + y;
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp," >>> LOCATION==%ld(fast-y, top-down)\n", location); }
        }
    } else {  // fast-X,  data
        if ( strcmp(uwpkfcvm_configuration->seek_axis, "fast-x") == 0 ||
                     strcmp(uwpkfcvm_configuration->seek_axis, "fast-X") == 0 ) { // fast-x,  uwpkfcvm 
            if(strcmp(uwpkfcvm_configuration->seek_direction, "bottom-up") == 0) { 
               location = ((long)z * uwpkfcvm_configuration->nx * uwpkfcvm_configuration->ny) + (y * uwpkfcvm_configuration->nx) + x;
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp," >>> LOCATION==%ld(fast-x, bottom-up)\n", location); }
                } else { // bottom-up
                    location = ((long)((uwpkfcvm_configuration->nz -1)- z) * uwpkfcvm_configuration->nx * uwpkfcvm_configuration->ny) + (y * uwpkfcvm_configuration->nx) + x;
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp," >>> LOCATION==%ld(fast-x, top-down)\n", location); }
            }
        }
    }

    // Check our loaded components of the model.
    if (uwpkfcvm_velocity_model->vp_status == 2) {
        // Read from memory.
        ptr = (float *)uwpkfcvm_velocity_model->vp;
        data->vp = ptr[location];
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp,"     FOUND : vp %f\n", data->vp); }
    } else if (uwpkfcvm_velocity_model->vp_status == 1) {
        // Read from file.
        fp = (FILE *)uwpkfcvm_velocity_model->vp;
        fseek(fp, location * sizeof(float), SEEK_SET);
        float temp;
        fread(&(temp), sizeof(float), 1, fp);
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp,"     FOUND : vp %f\n", temp); }
        data->vp=temp;
    }
    if (uwpkfcvm_velocity_model->vs_status == 2) {
        // Read from memory.
        ptr = (float *)uwpkfcvm_velocity_model->vs;
        data->vs = ptr[location];

if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp,"     FOUND : vs %f\n", data->vs); }

    } else if (uwpkfcvm_velocity_model->vs_status == 1) {
        // Read from file.
        fp = (FILE *)uwpkfcvm_velocity_model->vs;
        fseek(fp, location * sizeof(float), SEEK_SET);
        float temp;
        fread(&(temp), sizeof(float), 1, fp);
if(uwpkfcvm_ucvm_debug) {fprintf(stderrfp,"     FOUND : vs %f\n", temp); }
        data->vs=temp;
    }

    /* Calculate density */
    if (data->vp > 0.0) { data->rho=uwpkfcvm_calculate_density(data->vp); }
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void uwpkfcvm_trilinear_interpolation(double x_percent, double y_percent, double z_percent,
                             uwpkfcvm_properties_t *eight_points, uwpkfcvm_properties_t *ret_properties) {
    uwpkfcvm_properties_t *temp_array = calloc(2, sizeof(uwpkfcvm_properties_t));
    uwpkfcvm_properties_t *four_points = eight_points;

    uwpkfcvm_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

    // Now advance the pointer four "cvms5_properties_t" spaces.
    four_points += 4;

    // Another interpolation.
    uwpkfcvm_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

    // Now linearly interpolate between the two.
    uwpkfcvm_linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

    free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void uwpkfcvm_bilinear_interpolation(double x_percent, double y_percent, uwpkfcvm_properties_t *four_points, uwpkfcvm_properties_t *ret_properties) {
    uwpkfcvm_properties_t *temp_array = calloc(2, sizeof(uwpkfcvm_properties_t));
    uwpkfcvm_linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
    uwpkfcvm_linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
    uwpkfcvm_linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);
    free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void uwpkfcvm_linear_interpolation(double percent, uwpkfcvm_properties_t *x0, uwpkfcvm_properties_t *x1, uwpkfcvm_properties_t *ret_properties) {
    ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
    ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
    ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
    ret_properties->qp  = (1 - percent) * x0->qp  + percent * x1->qp;
    ret_properties->qs  = (1 - percent) * x0->qs  + percent * x1->qs;
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

    if (uwpkfcvm_velocity_model) free(uwpkfcvm_velocity_model);
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
            if (strcmp(key, "depth_interval") == 0) config->depth_interval = atof(value);
            if (strcmp(key, "seek_axis") == 0) sprintf(config->seek_axis, "%s", value);
            if (strcmp(key, "seek_direction") == 0) sprintf(config->seek_direction, "%s", value);
            if (strcmp(key, "interpolation") == 0) { 
                config->interpolation=0;
                if (strcmp(value,"on") == 0) config->interpolation=1;
            }
        }
    }

    // Have we set up all uwpkfcvm_configuration parameters?
    if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' || 
        config->seek_direction[0] == '\0' || config->seek_axis[0] == '\0' ||
        config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
        config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
        config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
        config->depth_interval == 0) {
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
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int uwpkfcvm_try_reading_model(uwpkfcvm_model_t *model) {
    double base_malloc = uwpkfcvm_configuration->nx * uwpkfcvm_configuration->ny * uwpkfcvm_configuration->nz * sizeof(float);
    int file_count = 0;
    int all_read_to_memory =0;
    char current_file[128];
    FILE *fp;

    // Let's see what data we actually have.
    sprintf(current_file, "%s/vp.dat", uwpkfcvm_data_directory);
    if (access(current_file, R_OK) == 0) {
        if( !too_big() ) { // only if fit
            model->vp = (float *)malloc(base_malloc);
            if (model->vp != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vp, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->vp_status = 2;
            } else {
              model->vp = fopen(current_file, "rb");
              model->vp_status = 1;
            }
        } else {
            model->vp = fopen(current_file, "rb");
            model->vp_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/vs.dat", uwpkfcvm_data_directory);
    if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
            model->vs = (float *) malloc(base_malloc);
        if (model->vs != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vs, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->vs_status = 2;
            } else {
              model->vs = fopen(current_file, "rb");
              model->vs_status = 1;
            }
        } else {
            model->vs = fopen(current_file, "rb");
            model->vs_status = 1;
                }
        file_count++;
    }

    if (file_count == 0)
        return FAIL;
    else if (file_count > 0 && all_read_to_memory != file_count)
        return SUCCESS;
    else
        return 2;
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
