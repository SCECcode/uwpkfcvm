/** 
  
  kdtree_util.h
  -- kdtree access code

**/

#ifndef KDTREE_UTIL_H
#define KDTREE_UTIL_H


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "proj.h"

#define EARTH_RADIUS_M 6371000.0

extern int uwpkfcvm_ucvm_debug;
extern int uwpkfcvm_ucvm_debug_detail;
extern FILE *stderrfp;

// reference lat/lon/depth/material properties
// for 'background expansion'
// pnts
typedef struct KDlld {
   double lat;
   double lon;
   double depth; 
   double vs;
   double vp;
} KDlld;

// base on KDlld 
// vpnts == in approximate utm
// in utm + depth 3D space
// where x is utm_e, y is utm_n, z is -depth
typedef struct KDVec3 {
    double x;
    double y;
    double z;
    int lldindex;
} KDVec3;

// nodes
typedef struct KDNode3 {
    KDVec3 *point;
    int axis;  // 0=x, 1=y, 2=z
    struct KDNode3 *left;
    struct KDNode3 *right;
} KDNode3;

typedef struct {
    double x;
    double y;
    double z;
    int lldindex; // original data index -> KDlld3
    int axis;
    int left;  // index of left child (-1 if none)
    int right; // index of right child (-1 if none)
} KDNode3Disk;


// this is creating convex hull in 2D
// using UTM 2D space
typedef struct KDVec2 {
    double utm_e;
    double utm_n;
    int lldindex;
} KDVec2;

/** access **/ 
void lld_to_xyz(KDVec3 *vp, double lat, double lon, double depth, int lldidx, PJ *_geo2utm);
void dump_v3pnts(KDVec3 *vp, int n);
KDNode3* build_v3kdtree(KDVec3 *pts, int n, int depth);
void free_v3kdtree(KDNode3 *node);
int flatten_v3kdtree(KDNode3 *node, KDNode3Disk *out, int *pos);
void write_flatten_v3kdtree(const char *fname, KDNode3Disk *nodes, int n);
KDNode3Disk *read_flatten_v3kdtree(const char *fname, int n);

void lld_to_en(KDVec2 *vp, KDlld *lld, int lldindex, PJ *_geo2utm);
void xyz_to_en(KDVec2 *vp, KDVec3 *xyz);
void dump_v2pnts(KDVec2 *vp, int n);
int create_boundary_hull(KDVec2 *pnts2, int n, KDVec2 **hull);

int setup_to_utm(PJ **_geo2utm, int MODEL_ZONE);
int to_utm(PJ *_geo2utm, double geo_lon, double geo_lat, double *utm_e, double *utm_n);
int to_geo(PJ *_geo2utm, double point_u, double point_v, double *lon, double *lat);

void find_xyz_latlon(KDlld *pnts, int lldindex, int nX, int nY);
void find_latlon(KDlld *pnts, int lldindex);

void lldindex_to_idx(int offset, int nx, int ny, int *xidx, int *yidx, int *zidx);
int idx_to_lldindex(int nx, int ny, int xidx, int yidx, int zidx);
KDVec3 *find_xyz_by_lldindex(KDVec3 *xyz, int n, int target);
double dist_sq(KDVec3* a, KDVec3* b);

/** usage **/
void kdtree_nearest(KDNode3 *node, KDVec3* query, KDVec3 **best, double *best_dist, int recursive);
int nearest_point(KDVec3 *points, int n, KDVec3 *query);
int point_in_convex(KDVec2 *poly, int n, KDVec2 p);

#endif // UWPKFCVM_H


