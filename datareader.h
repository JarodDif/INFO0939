#ifndef DATA_READER_H
#define DATA_READER_H

#include <stdio.h>

#include "common.h"

typedef struct grid_descriptor {
  unsigned int numnodesx;
  unsigned int numnodesy;
  unsigned int numnodesz;

  double xmin;
  double ymin;
  double zmin;

  double xmax;
  double ymax;
  double zmax;

} grid_descriptor_t;

typedef struct data {
  grid_descriptor_t grid;

  double *vals;

} data_t;


data_t *allocate_data(grid_descriptor_t grid);

data_t *read_data(FILE *fp, grid_descriptor_t grid,
                            unsigned int *idx, double *time);

FILE* open_datafile(grid_descriptor_t *grid, int *numsteps, char *filename);

void closest_index(data_t *data, double x, double y, double z, 
                                 int *mc, int *nc, int *pc);

#endif