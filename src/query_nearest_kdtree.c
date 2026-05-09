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
int uwpkfcvm_ucvm_debug=0;
int uwpkfcvm_ucvm_debug_detail=0;
FILE *stderrfp=NULL;
int is_rigid=1;

int NX=18;
int NY=15;
int NZ=9;

PJ *_geo2utm = NULL;
int MODEL_ZONE=10;

int main(int argc, char **argv)
{

  char datafile[KD_MAX_LINE];
  char queryfile[KD_MAX_LINE];
  char line[KD_MAX_LINE];
  int maxnum;
  FILE *fp;

  double lat, lon, depth, vs, vp;
  KDlld *pnts=NULL;
  int pnts_size=0;
  KDVec3 *v3pnts=NULL;  // e,n,depth
  KDNode3 *v3nodes=NULL; // tree v3nodes for vec3
  int *pnts_zero_depth;  // track which v3pnts has z of 0
                         
                    
  KDVec2 *v2pnts=NULL; // e,n
  KDVec2 *v2pnts_boundary=NULL; 

  int hull_size=0;
  KDVec2  *v2hull=NULL;

  if(argc != 4) { usage(); }
  if(sscanf(argv[1], "%d", &maxnum) != 1) { return 1; }
  strcpy(datafile, argv[2]);
  strcpy(queryfile, argv[3]);

  setup_to_utm(&_geo2utm, MODEL_ZONE);

  if(uwpkfcvm_ucvm_debug || uwpkfcvm_ucvm_debug_detail) {
      stderrfp = fopen("query_debug.log", "w+");
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
         pnts[numread].depth=depth*1000;
         pnts[numread].vs=vs;
         pnts[numread].vp=vp;
         // fillin KDVec3
         lld_to_xyz(&v3pnts[numread], lat, lon, depth*1000, numread, _geo2utm);   
         
         pnts_zero_depth[numread]=0;
         if(depth == 0 && is_rigid) {
             pnts_zero_depth[numread]=1;
             zero_depth_cnt++;
         }
         numread++;
      }
  }
  pnts_size=numread;
	  
  fclose(fp);

  if(is_rigid) {
      int r_idx=0;
      /** whole layer **/
      v2pnts = malloc( zero_depth_cnt * sizeof(KDVec2));
      for(int i=0; i< numread; i++) {
         if(pnts_zero_depth[i]) {
           int lldindex=v3pnts[i].lldindex;
           KDlld *lld = &pnts[lldindex];
           xyz_to_en(&v2pnts[r_idx],&v3pnts[i]);
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
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"BOUNDARY e(%lf) n(%lf)\n", v2pnts_boundary[k].utm_e, v2pnts_boundary[k].utm_n); }
      }
      
      hull_size=create_boundary_hull(v2pnts_boundary, boundary_cnt, &v2hull); 
      fprintf(stderr, "Convex Hull (%d points)\n", hull_size);
  }

  // kdtree for all the sorted grid points
  if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"==== kdtree with -- %d grid points and sorted v3pnts\n",numread); }
  v3nodes = build_v3kdtree(v3pnts, numread, 0);
  dump_v3pnts(v3pnts, numread);

  // QUERY against kdtree
  fp=fopen(queryfile,"r");
  int numquery=0;
  KDVec3 query_xyz;
  KDVec2 query_eu;
  KDlld query_lld;
  while (fgets(line, KD_MAX_LINE, fp) != NULL ) {
      if(line[0]=='#') continue;  // a comment line
      if(sscanf(line,"%lf %lf %lf", &lon, &lat, &depth) == 3) {
         KDVec3 *best;
         double best_dist = -1;
fprintf(stderr,"\n>>%lf %lf %lf \n", lon,lat,depth);
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\nQUERY point -> lat(%lf) lon(%lf) depth(%lf) \n", lat, lon, depth); }

         lld_to_xyz(&query_xyz, lat, lon, depth, 0/* don't care */, _geo2utm);

         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\nSEARCH with kdtree \n"); }
fprintf(stderr,"nearest point: kdtree SEARCH\n"); 
         kdtree_nearest(v3nodes, &query_xyz, &best, &best_dist, 1);
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"     >>>main: kdtree_nearest, best lldindex(%d), best dist(%lf)\n", best->lldindex, best_dist); }
         find_xyz_latlon(pnts, best->lldindex,NX,NY);
         fprintf(stderr,"  %lf %lf %lf\n",
                   pnts[best->lldindex].lon, pnts[best->lldindex].lat, pnts[best->lldindex].depth);

         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\nSEARCH with brute force \n"); }
fprintf(stderr,"nearest point: brute force  SEARCH\n"); 
         int rc=nearest_point(v3pnts, numread, &query_xyz);
         if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"     >>>main: nearest_point,  best lldindex(%d)\n", rc); }
         fprintf(stderr,"  %lf %lf %lf\n", pnts[rc].lon, pnts[rc].lat, pnts[rc].depth);

         if (is_rigid) {
fprintf(stderr,"in model: hull SEARCH\n"); 
           query_lld.lat=lat;
           query_lld.lon=lon;
           query_lld.depth=depth;
           lld_to_en(&query_eu, &query_lld, 0/* don't care */, _geo2utm);
           int yes_in=point_in_convex(v2hull, hull_size, query_eu);
           if(yes_in) {
              fprintf(stderr,"  in model\n");
	   } else { fprintf(stderr,"  out model\n");}
         }

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
  if(v2hull) free(v2hull);
  free(_geo2utm);
  if(uwpkfcvm_ucvm_debug) { fprintf(stderrfp,"\n..DONE..\n"); }
}
