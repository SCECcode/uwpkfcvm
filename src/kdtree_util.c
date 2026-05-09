/** 
  
  kdtree_util.c

**/

#include <math.h>
#include "kdtree_util.h"

void lld_to_xyz(KDVec3 *p, double lat, double lon, double depth, int lldidx, PJ *_geo2utm) {
    double utm_e;
    double utm_n;

    to_utm(_geo2utm, lon, lat, &utm_e, &utm_n);

    p->x = utm_e;
    p->y = utm_n;
    p->z = -depth;
    p->lldindex=lldidx;  // index into the whole data stream
    if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"    lld_to_xyz (%d), %lf %lf %lf(negative)\n", lldidx, p->x, p->y, p->z); }
}

// to ECEF global space -- global 3D Cartesian space
void lld_to_xyzECEF(KDVec3 *p, double lat, double lon, double depth, int lldidx)
{
    const double a = 6378137.0;            // WGS84 semi-major axis
    const double e2 = 6.69437999014e-3;    // eccentricity^2

    double nlat = lat * M_PI / 180.0;
    double nlon = lon * M_PI / 180.0;
    double h = depth * -1.0;  // depth → positive depth

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
       if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp," DUMP v3pnts : i(%d)  is lldindex(%d) \n",i, x->lldindex); }
     }
}

void dump_v2pnts(KDVec2 *vp, int n) {
     for(int i=0; i < n; i++) {
       KDVec2 *p = &vp[i];
       if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp," DUMP v2pnts : i(%d)  is lldindex(%d) \n",i, p->lldindex); }
     }
}

void _calc_idx(int lldindex, int *xidx, int *yidx, int *zidx, int nX, int nY) {
  int tmp= lldindex % (nX * nY);

  *zidx= lldindex / (nX * nY);
  *yidx= tmp / nX;
  *xidx= tmp % nX;
}

void find_xyz_latlon(KDlld *pnts, int lldindex,int nX, int nY) {
     // from lldindex figure out what is xidx,yidx,zidx
     int xidx; 
     int yidx; 
     int zidx; 
     _calc_idx(lldindex, &xidx, &yidx, &zidx, nX, nY);
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
    if(a->lldindex == b->lldindex)
      return 0.0;
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    double dz = a->z - b->z;
    double rc = sqrt(dx*dx + dy*dy + dz*dz);
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
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"   nearest_point : found (%d)\n", rc); }
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

    double d = dist_sq(node->point, query);
    if ( *best_dist == -1 || d < *best_dist) {
        *best_dist = d;
        *best = node->point;
	if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp," step kdtree_nearest(%d)  -- assign to best (%d) %lf\n", recursive, (*best)->lldindex, *best_dist); }
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
int setup_to_utm(PJ **_geo2utm, int MODEL_ZONE) {
  char _projstr[128];

  snprintf(_projstr, sizeof(_projstr), "+proj=utm +zone=%d +datum=NAD27 +units=m +no_defs", MODEL_ZONE);
  *_geo2utm=proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", _projstr, NULL);

  if(! *_geo2utm) {
    fprintf(stderr,"XX %s\n",(char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
    return 1;
  }

  *_geo2utm = proj_normalize_for_visualization(PJ_DEFAULT_CTX, *_geo2utm);

  return 0;
}

int to_utm(PJ *_geo2utm, double geo_lon, double geo_lat, double *utm_e, double *utm_n) {
  PJ_COORD xyzSrc = proj_coord(geo_lon, geo_lat, 0.0, HUGE_VAL);
  PJ_COORD xyzDest = proj_trans(_geo2utm, PJ_FWD, xyzSrc);
  int err = proj_context_errno(PJ_DEFAULT_CTX);
  if (err) {
     fprintf(stderr, "Error occurred while transforming latitude=%.4f, longitude=%.4f to UTM.\n", geo_lat, geo_lon);
     fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
     return 1;
  }
  *utm_e = xyzDest.xyzt.x;
  *utm_n = xyzDest.xyzt.y;
  return err;
}

int to_geo(PJ *_geo2utm, double point_u, double point_v, double *lon, double *lat) {
    PJ_COORD xyzSrc;
    xyzSrc.xyzt.x=point_u;
    xyzSrc.xyzt.y=point_v;
    PJ_COORD xyzDest = proj_trans(_geo2utm, PJ_INV, xyzSrc);

    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming u=%.4f, v=%.4f to Geo.\n",
              point_u, point_v);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return PJ_LOG_ERROR;
    }
    *lon=xyzDest.lp.lam;
    *lat=xyzDest.lp.phi;
    return err;
}


/*
 * Convert lld to ENU 
 */
void lld_to_en(KDVec2 *enu, KDlld *lld, int lldindex, PJ *_geo2utm) {
    double utm_e;
    double utm_n;

    to_utm(_geo2utm, lld->lon, lld->lat, &utm_e, &utm_n);

    enu->utm_e = utm_e;
    enu->utm_n = utm_n;

    enu->lldindex = lldindex;  // index into the whole data stream
}

void xyz_to_en(KDVec2 *enu, KDVec3 *xyz) {
    enu->utm_e = xyz->x;
    enu->utm_n = xyz->y;
    enu->lldindex = xyz->lldindex;  // index into the whole data stream
}


int cmp_vec2(const void *a, const void *b) {
    KDVec2 *p = (KDVec2 *)a;
    KDVec2 *q = (KDVec2 *)b;

    if (p->utm_e < q->utm_e) return -1;
    if (p->utm_e > q->utm_e) return 1;
    if (p->utm_n < q->utm_n) return -1;
    if (p->utm_n > q->utm_n) return 1;
    return 0;
}

double cross(KDVec2 O, KDVec2 A, KDVec2 B) {
    return (A.utm_e - O.utm_e)*(B.utm_n - O.utm_n) - (A.utm_n - O.utm_n)*(B.utm_e - O.utm_e);
}

int convex_hull(KDVec2 *pts, int n, KDVec2 *H) {
    if (n < 3) return 0;

    qsort(pts, n, sizeof(KDVec2), cmp_vec2);

    int k = 0;

    // lower hull
    for (int i = 0; i < n; i++) {
        while (k >= 2 && cross(H[k-2], H[k-1], pts[i]) <= 0) k--;
        H[k++] = pts[i];
    }

    // upper hull
    for (int i = n-2, t = k+1; i >= 0; i--) {
        while (k >= t && cross(H[k-2], H[k-1], pts[i]) <= 0) k--;
        H[k++] = pts[i];
    }

    return k-1;
}


// n is number of pts
int create_boundary_hull(KDVec2 *pnts2, int n, KDVec2 **hull) {

    *hull= (KDVec2 * ) malloc( (2*n) * sizeof (KDVec2) );

    int hsize = convex_hull(pnts2, n, *hull);

    if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp, "Convex Hull (%d points)\n", hsize); }
    for (int i = 0; i < hsize; i++) {
      if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"HULL: en(%f, %f)\n", (*hull)[i].utm_e, (*hull)[i].utm_n); }
    }
    
    return hsize;
}

/* ---------- Point-in-Convex ---------- */

int point_in_convex(KDVec2 *poly, int n, KDVec2 p) {
    if (cross(poly[0], poly[1], p) < 0) return 0;
    if (cross(poly[0], poly[n-1], p) > 0) return 0;

    int l = 1, r = n-1;

    while (r - l > 1) {
        int m = (l + r)/2;
        if (cross(poly[0], poly[m], p) >= 0)
            l = m;
        else
            r = m;
    }

    return cross(poly[l], poly[l+1], p) >= 0;
}

/**** interp with neighboring points *************/
// offset=(z_idx)*(ny * nx)+(y_idx)*(nx)+x_idx;
void lldindex_to_idx(int offset, int nx, int ny, int *xidx, int *yidx, int *zidx) {
  *zidx = offset / (nx * ny);
  int rem = offset % (nx * ny);

  *yidx = rem / nx;
  *xidx = rem % nx;
}

int idx_to_lldindex(int nx, int ny, int xidx, int yidx, int zidx) {
  int offset =(zidx)*(ny * nx)+(yidx)*(nx)+xidx;
  return offset;
}


// brute force look up
KDVec3 *find_xyz_by_lldindex(KDVec3 *xyz, int n, int target){

  KDVec3 *item;
  for(int i=0; i<n; i++) {
    item = &xyz[i];
    if(item->lldindex == target)
        return item;
  }
  return NULL;
}
