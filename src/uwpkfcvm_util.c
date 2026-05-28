/**

  uwpkfcvm_util.c

**/

#include <stdbool.h>

#include "uwpkfcvm.h"
#include "kdtree_util.h"

#define KD_MAX_LINE 1000

extern int uwpkfcvm_ucvm_debug;
extern int uwpkfcvm_ucvm_debug_detail;
extern FILE *stderrfp;
extern PJ *uwpkfcvm_geo2utm;

void setup_model(uwpkfcvm_model_t *model, int cnt) {
    // initialize all
    model->pnts = malloc(cnt * sizeof(KDlld));
    model->pnts_zero_depth = malloc(cnt * sizeof(int));
    model->zero_depth_cnt=0;

    model->v3pnts = malloc(cnt * sizeof(KDVec3));

    model->v2pnts = NULL;  // set later
    model->v2pnts_boundary = NULL;

    model->v2hull = NULL;
    model->v2hull_size = 0;

    model->tets = NULL;
    model->tets_cnt = 0;
}

void free_model(uwpkfcvm_model_t *model) {
    free(model->pnts);
    free(model->v3pnts);
    free(model->v2pnts);
    free(model->v2pnts_boundary);
    free(model->v2hull);   
    free(model->pnts_zero_depth);
    free(model->tets);
}

void load_model(uwpkfcvm_model_t *model, int NX, int NY, int NZ, FILE *fp) {
    int numread=0;
    char line[KD_MAX_LINE];
    int sz=NX * NY * NZ;

    setup_model(model, sz);

    double lat, lon, depth, vs, vp;

// load model data from external data file
    while (numread != sz && fgets(line, KD_MAX_LINE, fp) != NULL ) {
      if(line[0]=='#') continue;  // a comment line
      if (sscanf(line,"%lf %lf %lf %lf %lf", &lon, &lat, &depth, &vs, &vp) == 5) {
        model->pnts[numread].lat=lat;
        model->pnts[numread].lon=lon;
        model->pnts[numread].depth=depth * 1000; // in m
        model->pnts[numread].vs=vs;
        model->pnts[numread].vp=vp;
        // fillin KDVec3
        lld_to_xyz(&model->v3pnts[numread], lat, lon, (depth * 1000), numread);
        model->pnts_zero_depth[numread]=0;
        if(depth == 0) {
          model->pnts_zero_depth[numread]=1;
          model->zero_depth_cnt++;
        }
        numread++;
      }
    }
    model->pnts_size=numread;
    model->nx=NX;
    model->ny=NY;
    model->nz=NZ;

// setup the supporting structure for searching 
    model->v2pnts = malloc( model->zero_depth_cnt * sizeof(KDVec2));

    int r_idx=0;
    for(int i=0; i< numread; i++) {
      if(model->pnts_zero_depth[i]) {
	int lldindex=model->v3pnts[i].lldindex;
        KDlld *lld = &model->pnts[lldindex];
        lld_to_en(&model->v2pnts[r_idx],lld, lldindex, uwpkfcvm_geo2utm);
        r_idx++;
      }
    }
  

// collect up boundary pnts in order
    int boundary_cnt= (2*NX) + (2*NY) - 4;
    model->v2pnts_boundary= malloc( boundary_cnt * sizeof(KDVec2));

    int b_idx=0;
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"boundary_cnt %d\n",boundary_cnt); }
    // bottom row
    for(int i=0; i<NX; i++) {
      int t=i;
      model->v2pnts_boundary[b_idx]= model->v2pnts[t];
      if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"1 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
      b_idx++;
    }
    // right from bottom up
    for(int j=2; j<NY; j++) {  
      int t=(j * NX)-1;
      model->v2pnts_boundary[b_idx] = model->v2pnts[t];
      if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"2 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
      b_idx++;
    }
    // top row in reverse
    for(int i=1; i<NX; i++) { 
      int t= (NY * NX) - i;
      model->v2pnts_boundary[b_idx] = model->v2pnts[t];
        if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"3 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
        b_idx++;
    }
    // left from top down 
    for(int i=1; i<NY; i++) { 
      int t= (NX * (NY - i));
      model->v2pnts_boundary[b_idx] =  model->v2pnts[t];
      if(uwpkfcvm_ucvm_debug_detail) { fprintf(stderrfp,"4 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
      b_idx++;
    }
    // should be same as boundary_cnt
    model->boundary_size=b_idx;
 
    if(uwpkfcvm_ucvm_debug) { 
      fprintf(stderrfp,"== v2pnts (whole layer at depth 0) %d\n",r_idx);

      if(uwpkfcvm_ucvm_debug_detail) { 
        fprintf(stderrfp,"v2pnts (whole layer at depth 0) %d\n",r_idx);
        dump_v2pnts(model->v2pnts, r_idx );
      }
    }
  
    if(uwpkfcvm_ucvm_debug) { 
      fprintf(stderrfp,"== v2pnts_boundary (just boundary at depth 0) %d\n",b_idx);
      if(uwpkfcvm_ucvm_debug_detail) { 
        fprintf(stderrfp,"dump v2pnts_boundary (at depth 0) %d\n",b_idx);
        for(int k=0; k< boundary_cnt; k++) {
          find_latlon(model->pnts, model->v2pnts_boundary[k].lldindex);
          fprintf(stderrfp,"BOUNDARY e(%lf) n(%lf)\n",
	  		   model->v2pnts_boundary[k].utm_e, model->v2pnts_boundary[k].utm_n);
	}   
      }
    } 

    model->v2hull_size=create_boundary_hull(model->v2pnts_boundary, boundary_cnt, &model->v2hull);
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp, "Convex Hull (%d points)\n", model->v2hull_size); }

    if(build_tetrahedra_from_grid(model->nx, model->ny, model->nz, &model->tets, &model->tets_cnt) != 0) {
       fprintf(stderr, "BAD: fail to build Tet\n");
    }

// build kdtree for nearest neighbor searches
/* conflicting with the Tet interpolation
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"==== kdtree with -- %d grid points and sorted v3pnts\n",numread); }
    model->v3nodes = build_v3kdtree(model->v3pnts, numread, 0);
*/
    if(uwpkfcvm_ucvm_debug_detail) { dump_v3pnts(model->v3pnts, numread); }

}


int in_model(uwpkfcvm_model_t *model, double lat, double lon, double depth) {
    KDVec2 query_eu;
    KDlld query_lld;
    query_lld.lat=lat;
    query_lld.lon=lon;
    query_lld.depth=depth;
    lld_to_en(&query_eu, &query_lld, -1/* don't care */, uwpkfcvm_geo2utm);
    int rc=point_in_convex(model->v2hull, model->v2hull_size, query_eu);

    if(uwpkfcvm_ucvm_debug_detail) { 
      if(rc) { fprintf(stderrfp,"  in model\n"); }
      else { fprintf(stderrfp,"  out model\n");}
    }

    return rc;
}

// return the llindex into the model
int nearest_neighbor(uwpkfcvm_model_t *model, double lat, double lon, double depth, int total) {
    KDVec3 query_xyz;
    KDVec3 *best;
    double best_dist = -1;
    int rc;

    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\n..nearest_neighbor to target => lon(%f) lat(%f) depth(%f)\n", lon,lat,depth); }

    lld_to_xyz(&query_xyz, lat, lon, depth, -1/* don't care */);

/*
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"SEARCH with kdtree_nearest \n"); }
    kdtree_nearest(model->v3nodes, &query_xyz, &best, &best_dist, 1);
    rc=best->lldindex;
    if(uwpkfcvm_ucvm_debug) { 
	fprintf(stderrfp,"     >>>main: kdtree_nearest, best lldindex(%d), best dist(%lf)\n", best->lldindex, best_dist);
        fprintf(stderrfp,"  %lf %lf %lf\n", model->pnts[rc].lon, model->pnts[rc].lat, model->pnts[rc].depth);
    }
*/

    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"SEARCH with brute force \n"); }
    rc=nearest_point(model->v3pnts, total, &query_xyz);
    if(uwpkfcvm_ucvm_debug) { 
	fprintf(stderrfp,"     >>>main: nearest_point,  best lldindex(%d)\n", rc);
        fprintf(stderrfp,"  %lf %lf %lf\n", model->pnts[rc].lon, model->pnts[rc].lat, model->pnts[rc].depth);
    }

    return rc;
}


double vs_by_offset(uwpkfcvm_model_t *model, int loc) {
    if(loc == -1) { fprintf(stderr,"BAD.. bad vs location %d\n", loc); }
    return model->pnts[loc].vs;
}

double vp_by_offset(uwpkfcvm_model_t *model, int loc) {
    if(loc == -1) { fprintf(stderr,"BAD.. bad vp location %d\n", loc); }
    return model->pnts[loc].vp;
}


/**************** tetrahedral *******************/
KDVec3 _vsub(KDVec3 a, KDVec3 b) {
    KDVec3 r = {a.x - b.x, a.y - b.y, a.z - b.z, -1};
    return r;
}

KDVec3 _vcross(KDVec3 a, KDVec3 b) {
    KDVec3 r = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
	-1
    };
    return r;
}

/*
 * Returns barycentric coordinates (l0,l1,l2,l3) of p with respect to tetra
 * (p0,p1,p2,p3). Returns false if the tetra is degenerate.
 */
static bool _barycentric_tet(KDVec3 p0, KDVec3 p1, KDVec3 p2, KDVec3 p3, KDVec3 p,
                            double *l0, double *l1, double *l2, double *l3)
{
    KDVec3 a = _vsub(p1, p0);
    KDVec3 b = _vsub(p2, p0);
    KDVec3 c = _vsub(p3, p0);
    KDVec3 ap = _vsub(p,  p0);

    KDVec3 tmp= _vcross(b, c);
    double denom = dist_sq(&a, &tmp);
    if (fabs(denom) < 1e-20) {
        return false;
    }

    tmp= _vcross(b, c);
    double w1 = dist_sq(&ap, &tmp) / denom;
    tmp = _vcross(ap,c);
    double w2 = dist_sq(&a, &tmp)/ denom;
    tmp=_vcross(b, ap);
    double w3 = dist_sq(&a,  &tmp) / denom;
    double w0 = 1.0 - w1 - w2 - w3;

    if (l0) *l0 = w0;
    if (l1) *l1 = w1;
    if (l2) *l2 = w2;
    if (l3) *l3 = w3;

    return true;
}

static bool _point_in_tet(KDVec3 p0, KDVec3 p1, KDVec3 p2, KDVec3 p3, KDVec3 p,
                         double eps,
                         double *l0, double *l1, double *l2, double *l3) {
    if (!_barycentric_tet(p0, p1, p2, p3, p, l0, l1, l2, l3)) {
        return false;
    }

    return (*l0 >= -eps && *l1 >= -eps && *l2 >= -eps && *l3 >= -eps);
}


static double _tet_interp_scalar(KDVec3 q,
                                KDVec3 p0, KDVec3 p1, KDVec3 p2, KDVec3 p3,
                                double v0, double v1, double v2, double v3,
                                bool *ok) {
    double l0, l1, l2, l3;
    if (!_point_in_tet(p0, p1, p2, p3, q, 1e-12, &l0, &l1, &l2, &l3)) {
        if (ok) *ok = false;
        return NAN;
    }

    if (ok) *ok = true;
    return l0 * v0 + l1 * v1 + l2 * v2 + l3 * v3;
}

/*
/* Try the 6 tetrahedra of a hexahedral cell */
double interp_cell_tets(KDVec3 q, KDVec3 p[8], double v[8], bool *ok) {
    static const int tets[6][4] = {
        {0, 1, 3, 7},
        {0, 3, 2, 7},
        {0, 2, 6, 7},
        {0, 6, 4, 7},
        {0, 4, 5, 7},
        {0, 5, 1, 7}
    };

    for (int t = 0; t < 6; ++t) {
        int i0 = tets[t][0];
        int i1 = tets[t][1];
        int i2 = tets[t][2];
        int i3 = tets[t][3];

        double out = _tet_interp_scalar(q,
                                       p[i0], p[i1], p[i2], p[i3],
                                       v[i0], v[i1], v[i2], v[i3],
                                       ok);
        if (ok && *ok) {
            return out;
        }
    }

    if (ok) *ok = false;
    return NAN;
}

/* Interpolate over tetrahedra.
 *
 * model = model data structure
 * data  = returning data 
 * pts   = node coordinates
 * npts  = number of nodes
 * tets  = tetrahedra connectivity
 * ntets = number of tetrahedra
 * p     = query point
 *
 * Returns -1 if the point is outside all tetrahedra.
 */

int interpolate_tetra_mesh(uwpkfcvm_model_t *model, uwpkfcvm_properties_t *data,
            const KDVec3 *pts, int npts, const KDTet *tets, int ntets, KDVec3 p) {

    (void)npts; /* kept in case you want to add validation */
    const double eps = 1e-12;
    double vp_final;
    double vs_final;

    for (int t = 0; t < ntets; ++t) {
if(t==0) fprintf(stderr,"look into tet... %d(%d)\n", t,ntets);
        int i0 = tets[t].v[0];
        int i1 = tets[t].v[1];
        int i2 = tets[t].v[2];
        int i3 = tets[t].v[3];
if(t<5)
{
fprintf(stderr,"   %d %d %d %d\n",i0, i1, i2, i3);
print_latlon_by_lldindex(model->pnts,i0);
print_latlon_by_lldindex(model->pnts,i1);
print_latlon_by_lldindex(model->pnts,i2);
print_latlon_by_lldindex(model->pnts,i3);
fprintf(stderr,"\n");
}

        KDVec3 p0 = pts[i0];
        KDVec3 p1 = pts[i1];
        KDVec3 p2 = pts[i2];
        KDVec3 p3 = pts[i3];

        double l0, l1, l2, l3;
        if (_point_in_tet(p0, p1, p2, p3, p, eps, &l0, &l1, &l2, &l3)) {
fprintf(stderr,"FOUND IT... %d\n", t);
            vp_final=l0 * vp_by_offset(model,i0)+l1 * vp_by_offset(model,i1)+
                       l2 * vp_by_offset(model,i2)+l3 * vp_by_offset(model,i3);
            vs_final=l0 * vs_by_offset(model,i0)+l1 * vs_by_offset(model,i1)+
                        l2 * vs_by_offset(model,i2)+l3 * vp_by_offset(model,i3);

            data->vs=vs_final;
            data->vp=vp_final;
            /* Calculate density */
            if (data->vp > 0.0) { data->rho=uwpkfcvm_calculate_density(data->vp); }
        }
    }

    return -1;
}

