/** 
  
  kdtree_util.h
  -- kdtree access code

**/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define EARTH_RADIUS_M 6371000.0

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


// this is for building triangle surface from tthree KDVec2 points
typedef struct KDVec2 {
    double x;
    double y;
    int lldindex;
} KDVec2;

typedef struct KDTriangle {
    KDVec2 *a;
    KDVec2 *b;
    KDVec2 *c;
    KDVec2 *centroid;
} KDTriangle;


/** access **/ 
void lld_to_xyz(KDVec3 *vp, double lat, double lon, double depth, int lldidx);
void dump_v3pnts(KDVec3 *vp, int n);
KDNode3* build_v3kdtree(KDVec3 *pts, int n, int depth);
void free_v3kdtree(KDNode3 *node);
int flatten_v3kdtree(KDNode3 *node, KDNode3Disk *out, int *pos);
void write_flatten_v3kdtree(const char *fname, KDNode3Disk *nodes, int n);
KDNode3Disk *read_flatten_v3kdtree(const char *fname, int n);

void xyz_to_xy(KDVec2 *vp, KDVec3 *vs);
void dump_v2pnts(KDVec2 *vp, int n);


void find_xyz_latlon(KDlld *pnts, int lldindex);
void find_latlon(KDlld *pnts, int lldindex);

/** usage **/
double dist_sq(KDVec3* a, KDVec3* b);
void kdtree_nearest(KDNode3 *node, KDVec3* query, KDVec3 **best, double *best_dist, int recursive);
int nearest_point(KDVec3 *points, int n, KDVec3 *query);
