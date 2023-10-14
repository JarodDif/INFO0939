#include <stdlib.h>

#include "datareader.h"

data_t *allocate_data(grid_descriptor_t grid) {
  size_t numnodes = NUMNODESTOT(grid);
  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return NULL;
  }

  data_t *data;
  if ((data = malloc(sizeof(data_t))) == NULL ||
      (data->vals = malloc(numnodes * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return NULL;
  }

  data->grid = grid;

  return data;
}

data_t *read_data(FILE *fp, grid_descriptor_t grid, unsigned int *idx,
                  double *time) {
  if (fp == NULL) {
    DEBUG_PRINT("Invalid NULL file pointer");
    return NULL;
  }

  double ltime;
  unsigned int lidx;

  size_t numnodes = NUMNODESTOT(grid);

  data_t *data;
  if ((data = allocate_data(grid)) == NULL) {
    DEBUG_PRINT("Failed to allocate data");
    return NULL;
  }

  int readok = (fread(&lidx, sizeof(unsigned int), 1, fp) != 1 ||
                fread(&ltime, sizeof(double), 1, fp) != 1 ||
                fread(data->vals, sizeof(double), numnodes, fp) != numnodes);

  if (readok != 0) {
    DEBUG_PRINT("Failed to read data");
    return NULL;
  }

  if (idx != NULL)
    *idx = lidx;
  if (time != NULL)
    *time = ltime;

  return data;
}

FILE *open_datafile(grid_descriptor_t *grid, int *numsteps, char *filename) {
  if (grid == NULL || filename == NULL) {
    DEBUG_PRINT("Invalid NULL grid or filename");
    return NULL;
  }

  FILE *fp;
  if ((fp = fopen(filename, "rb")) == NULL) {
    DEBUG_PRINTF("Failed to open file '%s'", filename);
    return NULL;
  }

  fseek(fp, 0, SEEK_END);
  size_t filesize = ftell(fp);
  rewind(fp);

  int readok = (fread(&grid->numnodesx, sizeof(unsigned int), 1, fp) != 1 ||
                fread(&grid->numnodesy, sizeof(unsigned int), 1, fp) != 1 ||
                fread(&grid->numnodesz, sizeof(unsigned int), 1, fp) != 1 ||
                fread(&grid->xmin, sizeof(double), 1, fp) != 1 ||
                fread(&grid->xmax, sizeof(double), 1, fp) != 1 ||
                fread(&grid->ymin, sizeof(double), 1, fp) != 1 ||
                fread(&grid->ymax, sizeof(double), 1, fp) != 1 ||
                fread(&grid->zmin, sizeof(double), 1, fp) != 1 ||
                fread(&grid->zmax, sizeof(double), 1, fp) != 1);

  if (readok != 0) {
    DEBUG_PRINTF("Failed to read header of file '%s'", filename);
    fclose(fp);
    return NULL;
  }

  size_t numnodestot = grid->numnodesx * grid->numnodesy * grid->numnodesz;

  size_t datasize = filesize - 6 * sizeof(double) - 3 * sizeof(int);
  size_t perstepsize =
      numnodestot * sizeof(double) + sizeof(double) + sizeof(unsigned int);

  if (datasize % perstepsize != 0) {
    DEBUG_PRINTF("Data size is inconsistent with number of nodes (%lu, %lu)",
                 datasize, perstepsize);

    fclose(fp);
    return NULL;
  }

  if (numsteps != NULL) {
    *numsteps = (datasize / perstepsize);
  }

  return fp;
}

void closest_index(data_t *data, double x, double y, double z, int *mc, int *nc,
                   int *pc) {
  int m = (int)((x - XMIN(data)) / (XMAX(data) - XMIN(data)) * NUMNODESX(data));
  int n = (int)((y - YMIN(data)) / (YMAX(data) - YMIN(data)) * NUMNODESY(data));
  int p = (int)((z - ZMIN(data)) / (ZMAX(data) - ZMIN(data)) * NUMNODESZ(data));

  *mc = (m < 0) ? 0 : (m > NUMNODESX(data) - 1) ? NUMNODESX(data) - 1 : m;
  *nc = (n < 0) ? 0 : (n > NUMNODESY(data) - 1) ? NUMNODESY(data) - 1 : n;
  *pc = (p < 0) ? 0 : (p > NUMNODESZ(data) - 1) ? NUMNODESZ(data) - 1 : p;
}