/** 
  
  build_kdtree.c
  -- turn a list of lat-lon-z into a kdtree, 
     flatten the kdtree, 
     and write the flatten kdtree into a binary file
**/

#include <stdio.h>
#include <string.h>
#include "kdtree_util.h"

/*
    build_kdtree num-of-points file-of-surface-latlonz flatten-file
    ./build_kdtree 4 kdtree_surface.in kdtree_surface_flat.bin
*/ 

void usage() {
  printf("Usage: build_kdtree num_points kdtree_surface.in kdtree_surface_flat.bin\n\n");
  exit(0);
}

int uwpkfcvm_ucvm_debug=1;
int uwpkfcvm_ucvm_debug_detail=0;
FILE *stderrfp=NULL;


#define KD_MAX_LINE 1000

int main(int argc, char **argv)
{
  char datafile[KD_MAX_LINE];
  char flatfile[KD_MAX_LINE];
  char line[KD_MAX_LINE];
  int maxnum;
  FILE *fp;
  double lat, lon, depth, vs, vp;
  KDlld *pnts;
  KDVec3 *vpnts;
  KDNode3 *nodes; // tree nodes
  KDNode3Disk *dnodes; // tree nodes for disk

  if(argc != 4) { usage(); }
  if(sscanf(argv[1], "%d", &maxnum) != 1) { return 1; }
  strcpy(datafile, argv[2]);
  strcpy(flatfile, argv[3]);

  if(uwpkfcvm_ucvm_debug) {
      stderrfp = fopen("uwpkfcvm_debug.log", "w+");
      fprintf(stderrfp," ===== START ===== \n");
  }


  pnts = malloc(maxnum * sizeof(KDlld));
  vpnts = malloc(maxnum * sizeof(KDVec3));
  dnodes = malloc(maxnum * sizeof(KDNode3Disk));

  fp=fopen(datafile,"r");
  int numread=0;
// read all the points
  while (numread != maxnum && fgets(line, KD_MAX_LINE, fp) != NULL ) {
      if(line[0]=='#') continue;  // a comment line
				  //
      if (sscanf(line,"%lf %lf %lf %lf %lf", &lon, &lat, &depth, &vs, &vp) == 5) {
         pnts[numread].lat=lat;
         pnts[numread].lon=lon;
         pnts[numread].depth=depth;
         pnts[numread].vs=vs;
         pnts[numread].vp=vp;
	 lld_to_xyz(&vpnts[numread], lat, lon, depth, numread);   
         numread++;
      }
  }
  fclose(fp);
  nodes = build_v3kdtree(vpnts, numread, 0);
  fprintf(stderr,"kdtree with -- %d surface points\n",numread);

  int pos = 0;
  int flatten_top= flatten_v3kdtree(nodes, dnodes, &pos);

  fprintf(stderr,"flattend -- %d surface points\n",pos);

  write_flatten_v3kdtree(flatfile, dnodes, numread);

  // free all malloc
  free_v3kdtree(nodes);
  free(pnts);
  free(vpnts);
  fprintf(stderr,"\n..DONE..\n");
}
