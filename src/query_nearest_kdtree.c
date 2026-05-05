/** 
  
  query_nearest_kdtree.c
  -- convert a set of lat-lon-z into a kdtree, 
     make nearest point query against the kdtree
**/

#include <stdio.h>
#include <string.h>
#include "kdtree_util.h"

/*
    query_kdtree num-of-points file-of-surface-latlonz query_point
    ./query_kdtree 4 kdtree_surface.in kdtree_query.in
*/ 

void usage() {
  printf("Usage: query_kdtree num_points datafile query_list\n\n");
  exit(0);
}

#define KD_MAX_LINE 1000
int uwpkfcvm_ucvm_debug=1;
int uwpkfcvm_ucvm_debug_detail=0;
FILE *stderrfp=NULL;
int is_rigid=1;

int NX=18;
int NY=15;
int NZ=9;

int main(int argc, char **argv)
{

  char datafile[KD_MAX_LINE];
  char queryfile[KD_MAX_LINE];
  char line[KD_MAX_LINE];
  int maxnum;
  FILE *fp;
  double lat, lon, depth, vs, vp;
  KDlld *pnts=NULL;
  KDVec3 *v3pnts=NULL;  // x,y,z
  KDNode3 *v3nodes=NULL; // tree v3nodes for vec3
  int *pnts_zero_depth;  // track which v3pnts has z of 0
		    
  KDVec2 *v2pnts=NULL; // x,y
  KDVec2 *v2pnts_boundary=NULL; 

  if(argc != 4) { usage(); }
  if(sscanf(argv[1], "%d", &maxnum) != 1) { return 1; }
  strcpy(datafile, argv[2]);
  strcpy(queryfile, argv[3]);

  if(uwpkfcvm_ucvm_debug || uwpkfcvm_ucvm_debug_detail) {
      stderrfp = fopen("uwpkfcvm_debug.log", "w+");
      fprintf(stderrfp," ===== START ===== \n");
  }

  pnts = malloc(maxnum * sizeof(KDlld));
  v3pnts = malloc(maxnum * sizeof(KDVec3));
  pnts_zero_depth = malloc(maxnum * sizeof(int));

  fp=fopen(datafile,"r");
  int numread=0;
// read all the points
  int zero_depth_cnt=0; // if rigid then pick up points whose depth = 0
  if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"== BUILDING KDTREE ==\n"); }
  while (numread != maxnum && fgets(line, KD_MAX_LINE, fp) != NULL ) {
      if(line[0]=='#') continue;  // a comment line
      if (sscanf(line,"%lf %lf %lf %lf %lf", &lon, &lat, &depth, &vs, &vp) == 5) {
         pnts[numread].lat=lat;
         pnts[numread].lon=lon;
         pnts[numread].depth=depth;
         pnts[numread].vs=vs;
         pnts[numread].vp=vp;
         // fillin KDVec3
	 lld_to_xyz(&v3pnts[numread], lat, lon, depth, numread);   
	 
	 pnts_zero_depth[numread]=0;
         if(depth == 0 && is_rigid) {
             pnts_zero_depth[numread]=1;
	     zero_depth_cnt++;
         }
         numread++;
      }
  }
  fclose(fp);

  if(is_rigid) {
      int r_idx=0;
      /** whole layer **/
      v2pnts = malloc( zero_depth_cnt * sizeof(KDVec2));
      for(int i=0; i< numread; i++) {
         if(pnts_zero_depth[i]) {
           xyz_to_xy(&v2pnts[r_idx],&v3pnts[i]);
	   r_idx++;
         }
      }

      // in order
      int boundary_cnt= (2*NX) + (2*NY) - 4;
      v2pnts_boundary= malloc( boundary_cnt * sizeof(KDVec2));
      /** just boundary **/
      int b_idx=0;
      if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"boundary_cnt %d\n",boundary_cnt); }
      // bottom row
      for(int i=0; i<NX; i++) {
         int t=i;
         v2pnts_boundary[b_idx]= v2pnts[t];
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"1 add v2pnts_boundary %d\n",v2pnts[t].lldindex); }
         b_idx++;
      }
      // right from bottom up
      for(int j=2; j<NY; j++) {  
         int t=(j * NX)-1;
	 v2pnts_boundary[b_idx] =  v2pnts[t];
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"2 add v2pnts_boundary %d\n",v2pnts[t].lldindex); }
         b_idx++;
      }
      // top row in reverse
      for(int i=1; i<NX; i++) { 
         int t= (NY * NX) - i;
	 v2pnts_boundary[b_idx] =  v2pnts[t];
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"3 add v2pnts_boundary %d\n",v2pnts[t].lldindex); }
         b_idx++;
      }
      //
      // left from top down 
      for(int i=1; i<NY; i++) { 
         int t= (NX * (NY - i));
	 v2pnts_boundary[b_idx] =  v2pnts[t];
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"4 add v2pnts_boundary %d\n",v2pnts[t].lldindex); }
         b_idx++;
      }

      if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"dump v2pnts  %d\n",r_idx); }
      dump_v2pnts(v2pnts, r_idx );

      for(int k=0; k< boundary_cnt; k++) { 
         find_latlon(pnts, v2pnts_boundary[k].lldindex);
      }
  }

  // kdtree for all the sorted grid points
  if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"==== kdtree with -- %d grid points and sorted v3pnts\n",numread); }
  v3nodes = build_v3kdtree(v3pnts, numread, 0);
  dump_v3pnts(v3pnts, numread);

  // QUERY against kdtree
  fp=fopen(queryfile,"r");
  int numquery=0;
  KDVec3 query;
  while (fgets(line, KD_MAX_LINE, fp) != NULL ) {
      if(line[0]=='#') continue;  // a comment line
      if(sscanf(line,"%lf %lf %lf", &lon, &lat, &depth) == 3) {
         KDVec3 *best;
         double best_dist = -1;
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\nQUERY point -> lat(%lf) lon(%lf) depth(%lf) \n", lat, lon, depth); }
	 lld_to_xyz(&query, lat, lon, depth, 0/* don't care */);

	 
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\nSEARCH with kdtree \n"); }
	 kdtree_nearest(v3nodes, &query, &best, &best_dist, 1);
	 if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"     >>>main: kdtree_nearest, best lldindex(%d), best dist(%lf)\n", best->lldindex, best_dist); }
	 find_xyz_latlon(pnts, best->lldindex);

         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\nSEARCH with brute force \n"); }
	 int rc=nearest_point(v3pnts, numread, &query);
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"     >>>main: nearest_point,  best lldindex(%d)\n", rc); }
         numquery++;
      }
  }
  if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"query for %d locations\n", numquery); }

  // remember to free all malloc
  free_v3kdtree(v3nodes);
  free(pnts);
  free(v3pnts);
  if(v2pnts) free(v2pnts);
  if(v2pnts_boundary) free(v2pnts_boundary);
  if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\n..DONE..\n"); }
}
