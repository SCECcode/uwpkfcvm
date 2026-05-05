/** 
  
  kdtree_util.c.c

**/

#include <math.h>
#include "kdtree_util.h"

extern int uwpkfcvm_ucvm_debug;
extern int uwpkfcvm_ucvm_debug_detail;
extern FILE *stderrfp;
//
// number of Xs and Ys and Zs in the 3D grid space
int nX=18;
int nY=15;
int nZ=9;


// to ECEF global space
void lld_to_xyz(KDVec3 *p, double lat, double lon, double depth, int lldidx)
{
    const double a = 6378137.0;            // WGS84 semi-major axis
    const double e2 = 6.69437999014e-3;    // eccentricity^2

    double nlat = lat * M_PI / 180.0;
    double nlon = lon * M_PI / 180.0;
    double h = depth * 1000.0;  // depth → positive height

    double sin_lat = sin(nlat);
    double cos_lat = cos(nlat);

    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);

    p->x = (N + h) * cos_lat * cos(nlon);
    p->y = (N + h) * cos_lat * sin(nlon);
    p->z = (N * (1 - e2) + h) * sin_lat;
    p->lldindex=lldidx;  // index into the whole data stream

    if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"    lld_to_xyz (%d), %lf %lf %lf\n", lldidx, p->x, p->y, p->z); }
}

void dump_v3pnts(KDVec3 *vp, int n) {
     for(int i=0; i < n; i++) {
       KDVec3 *x = &vp[i];
       if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp," DUMP v3pnts : %d  is lldindex(%d) \n",i, x->lldindex); }
     }
}

void dump_v2pnts(KDVec2 *vp, int n) {
     for(int i=0; i < n; i++) {
       KDVec2 *p = &vp[i];
       if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp," DUMP v2pnts : %d  is lldindex(%d) \n",i, p->lldindex); }
     }
}

void _calc_idx(int lldindex, int *xidx, int *yidx, int *zidx) {
  int tmp= lldindex % (nX * nY);

  *zidx= lldindex / (nX * nY);
  *yidx= tmp / nX;
  *xidx= tmp % nX;
}

void find_xyz_latlon(KDlld *pnts, int lldindex) {
     // from lldindex figure out what is xidx,yidx,zidx
     int xidx; 
     int yidx; 
     int zidx; 
     _calc_idx(lldindex, &xidx, &yidx, &zidx);
     if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   FOUND_latlon: lldindex(%d), lon(%lf) lat(%lf) depth(%lf)\n              vs(%lf) vp(%lf)\n              xidx(%d) yidx(%d) zidx(%d)\n", 
		     lldindex, pnts[lldindex].lon, pnts[lldindex].lat, pnts[lldindex].depth,
		     pnts[lldindex].vs, pnts[lldindex].vp,
		     xidx, yidx, zidx); }
}
void find_latlon(KDlld *pnts, int lldindex) {
     if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"%lf %lf %lf\n", 
		      pnts[lldindex].lon, pnts[lldindex].lat, pnts[lldindex].depth);}
}

int cmp_x(const void *a, const void *b) { return ((KDVec3*)a)->x > ((KDVec3*)b)->x; }
int cmp_y(const void *a, const void *b) { return ((KDVec3*)a)->y > ((KDVec3*)b)->y; }
int cmp_z(const void *a, const void *b) { return ((KDVec3*)a)->z > ((KDVec3*)b)->z; }

// n=number of points
// depth start at 0
KDNode3* build_v3kdtree(KDVec3 *pts, int n, int depth)
{
    if (n <= 0) return NULL;

    int axis = depth % 3;
    if (axis == 0) qsort(pts, n, sizeof(KDVec3), cmp_x);
    if (axis == 1) qsort(pts, n, sizeof(KDVec3), cmp_y);
    if (axis == 2) qsort(pts, n, sizeof(KDVec3), cmp_z);

    int mid = n / 2;

    KDNode3 *node = malloc(sizeof(KDNode3));
    node->point = &pts[mid];
    node->axis = axis;

    node->left  = build_v3kdtree(pts, mid, depth+1);
    node->right = build_v3kdtree(pts+mid+1, n-mid-1, depth+1);

    return node;
}

void free_v3kdtree(KDNode3 *node) {
    if(node->left) free_v3kdtree(node->left);
    if(node->right) free_v3kdtree(node->right);
    free(node);
}

int flatten_v3kdtree(KDNode3 *node, KDNode3Disk *out, int *pos) {
    if (!node) return -1;

    int my_id = (*pos)++;

    out[my_id].x = node->point->x;
    out[my_id].y = node->point->y;
    out[my_id].z = node->point->z;
    out[my_id].axis = node->axis;
    out[my_id].lldindex = node->point->lldindex;

    out[my_id].left  = flatten_v3kdtree(node->left,  out, pos);
    out[my_id].right = flatten_v3kdtree(node->right, out, pos);

    return my_id;
}

void write_flatten_v3kdtree(const char *fname, KDNode3Disk *nodes, int n) {

    FILE *f = fopen(fname, "wb");
    if (!f) return;

    fwrite(&n, sizeof(int), 1, f);
    fwrite(nodes, sizeof(KDNode3Disk), n, f);

    fclose(f);
}

// read back kdtree from binary file
KDNode3Disk *read_flatten_v3kdtree(const char *fname, int n) {

    FILE *fp = fopen(fname, "rb");
    if (!fp) return NULL;

    KDNode3Disk *dnodes = malloc(n * sizeof(KDNode3Disk));
    fread(dnodes, sizeof(KDNode3Disk), n, fp);
    fclose(fp);

    return dnodes;
}

// XXX  rebuild structure ??


// nearest code
double dist_sq(KDVec3* a, KDVec3* b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;
    double rc = dx*dx + dy*dy + dz*dz;
    if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"              dist_sq : distance to (%d) is  a=%lf\n", a->lldindex, rc); }
     
    return rc;
}

int nearest_point(KDVec3 *points, int n, KDVec3 *query) {

    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"    ==CALLING nearest_point ==\n"); }
    int best_idx = -1;
    double best_dist = 1e100;

    for (int i = 0; i < n; i++) {
        double d = dist_sq(&points[i], query);
        if (d < best_dist) {
            best_dist = d;
            best_idx = i;
        }
    }
    int rc = points[best_idx].lldindex;
    return rc;
}


//t3D query = lld_to_xyz(q_lat, q_lon, q_depth, -1);
//Point3D best;
//double best_dist = DBL_MAX;
//kdtree_nearest(root, query, &best, &best_dist);
//printf("Nearest index: %d\n", best.index);
//This is recursive
void kdtree_nearest(KDNode3 *node, KDVec3 *query, KDVec3 **best, double *best_dist, int recursive) {
    if (!node) return;

    if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"    ==CALLING kdtree_nearest==\n"); }
    double d = dist_sq(node->point, query);
    if ( *best_dist == -1 || d < *best_dist) {
        *best_dist = d;
        *best = node->point;
	if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   == CALLING kdtree_nearest(%d)  -- assign to best (%d) %lf\n", recursive, (*best)->lldindex, *best_dist); }
    }

    int axis = node->axis;
    double diff;

    if (axis == 0) diff = query->x - node->point->x;
    else if (axis == 1) diff = query->y - node->point->y;
    else diff = query->z - node->point->z;

    KDNode3 *first  = diff < 0 ? node->left  : node->right;
    KDNode3 *second = diff < 0 ? node->right : node->left;

    kdtree_nearest(first, query, best, best_dist, recursive+1);

    if (diff*diff < *best_dist) {
        kdtree_nearest(second, query, best, best_dist, recursive+1);
    }
}



/************************************************************/
void xyz_to_xy(KDVec2 *p, KDVec3 *s) {
    p->x = s->x;
    p->y = s->y;
    p->lldindex = s->lldindex;  // index into the whole data stream
}

