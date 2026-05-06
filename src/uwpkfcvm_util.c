/**

  uwpkfcvm_util.c

**/

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
}

void free_model(uwpkfcvm_model_t *model) {
    free(model->pnts);
    free(model->v3pnts);
    free(model->v2pnts);
    free(model->v2pnts_boundary);
    free(model->v2hull);   
    free(model->pnts_zero_depth);
}

int load_model(uwpkfcvm_model_t *model, int NX, int NY, int NZ, FILE *fp) {
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
        model->pnts[numread].depth=depth;
        model->pnts[numread].vs=vs;
        model->pnts[numread].vp=vp;
        // fillin KDVec3
        lld_to_xyz(&model->v3pnts[numread], lat, lon, depth, numread);
  
        model->pnts_zero_depth[numread]=0;
        if(depth == 0) {
          model->pnts_zero_depth[numread]=1;
          model->zero_depth_cnt++;
        }
        numread++;
      }
    }

// setup the supporting structure for searching 
    model->v2pnts = malloc( model->zero_depth_cnt * sizeof(KDVec2));

    int r_idx=0;
    for(int i=0; i< numread; i++) {
      if(model->pnts_zero_depth[i]) {
        int lldindex=model->v3pnts[i].lldindex;
        KDlld *lld = &model->pnts[lldindex];
        lld_to_en(&model->v2pnts[r_idx],lld,lldindex,uwpkfcvm_geo2utm);
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
      if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"1 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
      b_idx++;
    }
    // right from bottom up
    for(int j=2; j<NY; j++) {  
      int t=(j * NX)-1;
      model->v2pnts_boundary[b_idx] = model->v2pnts[t];
      if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"2 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
      b_idx++;
    }
    // top row in reverse
    for(int i=1; i<NX; i++) { 
      int t= (NY * NX) - i;
      model->v2pnts_boundary[b_idx] = model->v2pnts[t];
        if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"3 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
        b_idx++;
    }
    // left from top down 
    for(int i=1; i<NY; i++) { 
      int t= (NX * (NY - i));
      model->v2pnts_boundary[b_idx] =  model->v2pnts[t];
      if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"4 add v2pnts_boundary %d\n",model->v2pnts[t].lldindex); }
      b_idx++;
    }
    // should be same as boundary_cnt
    model->boundary_size=b_idx;
 
    if(uwpkfcvm_ucvm_debug) { 
      fprintf(stderrfp,"dump v2pnts  %d\n",r_idx);
      dump_v2pnts(model->v2pnts, r_idx );
    }
  
    for(int k=0; k< boundary_cnt; k++) {
      if(uwpkfcvm_ucvm_debug) { find_latlon(model->pnts, model->v2pnts_boundary[k].lldindex); }
      if(uwpkfcvm_ucvm_debug) { 
        find_latlon(model->pnts, model->v2pnts_boundary[k].lldindex);
        if(uwpkfcvm_ucvm_debug) { 
          fprintf(stderrfp,"BOUNDARY e(%lf) n(%lf)\n",
			   model->v2pnts_boundary[k].utm_e, model->v2pnts_boundary[k].utm_n);
	}
      }
 
      model->v2hull_size=create_boundary_hull(model->v2pnts_boundary, boundary_cnt, &model->v2hull);
      if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp, "Convex Hull (%d points)\n", model->v2hull_size); }
    }

// build kdtree for nearest neighbor searches
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"==== kdtree with -- %d grid points and sorted v3pnts\n",numread); }
    model->v3nodes = build_v3kdtree(model->v3pnts, numread, 0);
    if(uwpkfcvm_ucvm_debug) { dump_v3pnts(model->v3pnts, numread); }
}


int in_model(uwpkfcvm_model_t *model, double lat, double lon, double depth) {
    KDVec2 query_eu;
    KDlld query_lld;
    query_lld.lat=lat;
    query_lld.lon=lon;
    query_lld.depth=depth;
    lld_to_en(&query_eu, &query_lld, 0/* don't care */, uwpkfcvm_geo2utm);
    int rc=point_in_convex(model->v2hull, model->v2hull_size, query_eu);

    if(uwpkfcvm_ucvm_debug) { 
      if(rc) { fprintf(stderrfp,"  in model\n"); }
      else { fprintf(stderrfp,"  out model\n");}
    }

    return rc;
}

// return the llindex into the model
int nearest_neighbor(uwpkfcvm_model_t *model, double lat, double lon, double depth) {
    KDVec3 query_xyz;
    KDVec3 *best;
    double best_dist = -1;

    lld_to_xyz(&query_xyz, lat, lon, depth, 0/* don't care */);
    kdtree_nearest(model->v3nodes, &query_xyz, &best, &best_dist, 1);
    if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"     >>>main: kdtree_nearest, best lldindex(%d), best dist(%lf)\n", best->lldindex, best_dist); }
    return best->lldindex;
}


double vs_by_location(uwpkfcvm_model_t *model, int loc) {
    return model->pnts[loc].vs;
}

double vp_by_location(uwpkfcvm_model_t *model, int loc) {
    return model->pnts[loc].vp;
}

