#include <mpi.h>

#include "fdtd.h"

// Global topology values that once set never change
MPI_Comm cart_comm;
int cart_rank;
int dims[3];
int coords[3];
int neighbors[6];

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  
  if (argc < 2) {
    printf("\nUsage: ./fdtd <param_file>\n\n");
    exit(1);
  }

  int world_size, rank, reorder = 0;
  int periods[3] = {0,0,0};

  // Set up topology
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Dims_create(world_size, 3, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &cart_comm);
  MPI_Comm_rank(cart_comm, &cart_rank);

  MPI_Cart_coords(cart_comm, cart_rank, 3, coords);
  MPI_Cart_shift(cart_comm, 0, 1, 
    &neighbors[LEFT], &neighbors[RIGHT]);
  MPI_Cart_shift(cart_comm, 1, 1, 
    &neighbors[FRONT], &neighbors[BACK]);
  MPI_Cart_shift(cart_comm, 2, 1, 
		&neighbors[DOWN], &neighbors[UP]);

  process_simulation_data_t simdata;
  init_simulation(&simdata, argv[1]);

#if 0
  int numtimesteps = floor(simdata.params.maxt / simdata.params.dt);

  double start = GET_TIME();
  for (int tstep = 0; tstep <= numtimesteps; tstep++) {
    apply_source(&simdata, tstep);

    if (simdata.params.outrate > 0 && (tstep % simdata.params.outrate) == 0) {
      for (int i = 0; i < simdata.params.numoutputs; i++) {
        data_t *output_data = NULL;

        switch (simdata.params.outputs[i].source) {
        case PRESSURE:
          output_data = simdata.pold;
          break;
        case VELOCITYX:
          output_data = simdata.vxold;
          break;
        case VELOCITYY:
          output_data = simdata.vyold;
          break;
        case VELOCITYZ:
          output_data = simdata.vzold;
          break;

        default:
          break;
        }

        double time = tstep * simdata.params.dt;
        write_output(&simdata.params.outputs[i], output_data, tstep, time);
      }
    }

    if (tstep > 0 && tstep % (numtimesteps / 10) == 0) {
      printf("step %8d/%d", tstep, numtimesteps);

      if (tstep != numtimesteps) {
        double elapsed_sofar = GET_TIME() - start;
        double timeperstep_sofar = elapsed_sofar / tstep;

        double eta = (numtimesteps - tstep) * timeperstep_sofar;

        printf(" (ETA: %8.3lf seconds)", eta);
      }

      printf("\n");
      fflush(stdout);
    }

    update_pressure(&simdata);
    update_velocities(&simdata);
    swap_timesteps(&simdata);
  }

  double elapsed = GET_TIME() - start;
  double numupdates =
      (double)NUMNODESTOT(simdata.pold->grid) * (numtimesteps + 1);
  double updatespers = numupdates / elapsed / 1e6;

  printf("\nElapsed %.6lf seconds (%.3lf Mupdates/s)\n\n", elapsed,
         updatespers);

#endif

  finalize_simulation(&simdata);

  MPI_Finalize();

  return 0;
}

/******************************************************************************
 * Utilities functions                                                        *
 ******************************************************************************/
char *copy_string(char *str) {
  size_t len;
  if (str == NULL || (len = strlen(str)) == 0) {
    DEBUG_PRINT("NULL of zero length string passed as argument");
    return NULL;
  }

  char *cpy;
  if ((cpy = malloc((len + 1) * sizeof(char))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return NULL;
  }

  return strcpy(cpy, str);
}

void closest_index(grid_t *grid, double x, double y, double z, int *cx, int *cy, int *cz) {
  int m = (int)((x - grid->xmin) / (grid->xmax - grid->xmin) * grid->numnodesx);
  int n = (int)((y - grid->ymin) / (grid->ymax - grid->ymin) * grid->numnodesy);
  int p = (int)((z - grid->zmin) / (grid->zmax - grid->zmin) * grid->numnodesz);

  *cx = (m < 0) ? 0 : (m > grid->numnodesx - 1) ? grid->numnodesx - 1 : m;
  *cy = (n < 0) ? 0 : (n > grid->numnodesy - 1) ? grid->numnodesy - 1 : n;
  *cz = (p < 0) ? 0 : (p > grid->numnodesz - 1) ? grid->numnodesz - 1 : p;
}

double trilinear_interpolation(data_t *data, double x, double y, double z) {
  //Compute all needed values
  double m = (x - XMIN(data)) / (XMAX(data) - XMIN(data)) * (NUMNODESX(data)-1);
  double n = (y - YMIN(data)) / (YMAX(data) - YMIN(data)) * (NUMNODESY(data)-1);
  double p = (z - ZMIN(data)) / (ZMAX(data) - ZMIN(data)) * (NUMNODESZ(data)-1);

  register int m0 = (int)m, n0 = (int)n, p0 = (int)p;
  m0 = (m0 < 0) ? 0 : (m0 > NUMNODESX(data) - 1) ? NUMNODESX(data) - 1 : m0;
  n0 = (n0 < 0) ? 0 : (n0 > NUMNODESY(data) - 1) ? NUMNODESY(data) - 1 : n0;
  p0 = (p0 < 0) ? 0 : (p0 > NUMNODESZ(data) - 1) ? NUMNODESZ(data) - 1 : p0;
  register double dm = m - m0, dn = n - n0, dp = p - p0;

  double c[8];

  for(int i=0; i < 8; i++){
    c[i] = GETVALUE(data, 
      (m0 != NUMNODESX(data)-1 && i&4)?(m0+1):m0,
      (n0 != NUMNODESY(data)-1 && i&2)?(n0+1):n0,
      (p0 != NUMNODESZ(data)-1 && i&1)?(p0+1):p0);
  }

  //reduce to cube to square
  for(int i=0; i < 4; i++){
    c[i] = c[i] * (1 - dm) + c[i+4] * dm;
  }

  //reduce square to line 
  for(int i=0; i < 2; i++){
    c[i] = c[i] * (1 - dn) + c[i+2] * dn;
  }

  //reduce line to point
  return c[0] * (1 - dp) + c[1] * dp;
}

void print_source(source_t *source) {
  printf(" Source infos:\n\n");

  if (source->type == AUDIO) {
    double duration = (double)source->numsamples / source->sampling;

    printf("          type: audio data file\n");
    printf("      sampling: %d Hz\n", source->sampling);
    printf("      duration: %g\n", duration);

  } else {
    printf("          type: sine wave\n");
    printf("     frequency: %g Hz\n", source->data[0]);
  }

  printf("    position x: %g\n", source->posx);
  printf("    position y: %g\n", source->posy);
  printf("    position z: %g\n\n", source->posy);
}

void print_output(output_t *output) {
  switch (output->source) {
  case PRESSURE:
    printf("      pressure: ");
    break;
  case VELOCITYX:
    printf("    velocity X: ");
    break;
  case VELOCITYY:
    printf("    velocity Y: ");
    break;
  case VELOCITYZ:
    printf("    velocity Z: ");
    break;

  default:
    break;
  }

  switch (output->type) {
  case ALL:
    printf("complete dump");
    break;
  case CUTX:
    printf("cut along the x axis at %g", output->posx);
    break;
  case CUTY:
    printf("cut along the y axis at %g", output->posy);
    break;
  case CUTZ:
    printf("cut along the z axis at %g", output->posz);
    break;
  case POINT:
    printf("single point at %g %g %g", output->posx, output->posy,
           output->posz);
    break;

  default:
    break;
  }

  printf(" to file %s\n", output->filename);
}

int process_setvalue(process_data_t *pdata, int m, int n, int p, double val) {
  //Verify we are inside the boundaries of the whole grid
  if(m < 0 || n < 0 || p < 0) 
    return 0.0;
  if(m >= pdata->grid.gnumx || n >= pdata->grid.gnumy || p >= pdata->grid.gnumz)
    return 0.0;

  int mbar = m - STARTM(pdata), nbar = n - STARTN(pdata), pbar = p - STARTP(pdata);
  int pnumnodesx = PNUMNODESX(pdata), pnumnodesy = PNUMNODESY(pdata), pnumnodesz = PNUMNODESZ(pdata);

  //Verify we are in the boundaries of the current subdomain
  if(mbar < -1 || mbar > pnumnodesx ||
      nbar < -1 || nbar > pnumnodesy ||
      pbar < -1 || pbar > pnumnodesz) {
    return 1;
  }

  if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
      nbar >= 0 && nbar <= pnumnodesy - 1 &&
      pbar >= 0 && pbar <= pnumnodesz - 1) {
    pdata->vals[pbar*pnumnodesx*pnumnodesy + nbar*pnumnodesx + mbar] = val;
  }
  else if(nbar >= 0 && nbar <= pnumnodesy - 1 &&
          pbar >= 0 && pbar <= pnumnodesz - 1) {
    pdata->ghostvals[(mbar == -1)?LEFT:RIGHT][pbar*PNUMNODESY(pdata) + nbar] = val;
  }
  else if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
          pbar >= 0 && pbar <= pnumnodesz - 1) {
    pdata->ghostvals[(nbar == -1)?FRONT:BACK][pbar*pnumnodesx + mbar] = val;
  }
  else if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
          nbar >= 0 && nbar <= pnumnodesy - 1) {
    pdata->ghostvals[(pbar == -1)?DOWN:UP][nbar*pnumnodesx + mbar] = val;
  }
  else{ return 1; }

  return 0;
}

double process_getvalue(process_data_t *pdata, int m, int n, int p) {
  //Verify we are inside the boundaries
  if(m < 0 || n < 0 || p < 0) 
    return 0.0;
  if(m >= pdata->grid.gnumx || n >= pdata->grid.gnumy || p >= pdata->grid.gnumz)
    return 0.0;
  
  int mbar = m - STARTM(pdata), nbar = n - STARTN(pdata), pbar = p - STARTP(pdata);
  int pnumnodesx = PNUMNODESX(pdata), pnumnodesy = PNUMNODESY(pdata), pnumnodesz = PNUMNODESZ(pdata);

  //Verify we are in the boundaries of the current subdomain
  if(mbar < -1 || mbar > pnumnodesx ||
      nbar < -1 || nbar > pnumnodesy ||
      pbar < -1 || pbar > pnumnodesz) {
    return 0.0;
  }

  if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
    nbar >= 0 && nbar <= pnumnodesy - 1 && 
    pbar >= 0 && pbar <= pnumnodesz - 1){
    return pdata->vals[pbar*pnumnodesx*pnumnodesy + nbar*pnumnodesx + mbar];
  }
  else if(nbar >= 0 && nbar <= pnumnodesy - 1 &&
          pbar >= 0 && pbar <= pnumnodesz - 1) {
    return pdata->ghostvals[(mbar == -1)?LEFT:RIGHT][pbar*PNUMNODESY(pdata) + nbar];
  }
  else if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
          pbar >= 0 && pbar <= pnumnodesz - 1) {
    return pdata->ghostvals[(nbar == -1)?FRONT:BACK][pbar*pnumnodesx + mbar];
  }
  else if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
          nbar >= 0 && nbar <= pnumnodesy - 1) {
    return pdata->ghostvals[(pbar == -1)?DOWN:UP][nbar*pnumnodesx + mbar];
  }

  return 0.0;
}

/******************************************************************************
 * Data functions                                                             *
 ******************************************************************************/
data_t* allocate_data(grid_t *grid){
  size_t numnodes = NUMNODESTOT(*grid);
  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return NULL;
  }

  data_t *data;
  if ((data = malloc(sizeof(data_t))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(data);
    return NULL;
  }

  if ((data->vals = malloc(numnodes * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(data->vals);
    free(data);
    return NULL;
  }

  data->grid = *grid;

  return data;
}

void fill_data(data_t *data, double value) {
  if (data == NULL) {
    DEBUG_PRINT("Invalid NULL data");
    return;
  }

  for (int p = 0; p < NUMNODESZ(data); p++) {
    for (int n = 0; n < NUMNODESY(data); n++) {
      for (int m = 0; m < NUMNODESX(data); m++) {
        SETVALUE(data, m, n, p, value);
      }
    }
  }
}

process_data_t *allocate_pdata(process_grid_t *grid) {
  int pnumnodesx = grid->endm - grid->startm + 1, 
    pnumnodesy = grid->endn - grid->startn + 1, 
    pnumnodesz = grid->endp - grid->startp + 1;
  size_t numnodes = pnumnodesx*pnumnodesy*pnumnodesz;
  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return NULL;
  }

  process_data_t *pdata;
  if ((pdata = malloc(sizeof(process_data_t))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(pdata);
    return NULL;
  }

  if ((pdata->vals = malloc(numnodes * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    free(pdata->vals);
    free(pdata);
    return NULL;
  }

  if((pdata->ghostvals = malloc(NEIGHBOR_TYPE_END * sizeof(double*))) == NULL){
    DEBUG_PRINT("Failed to allocate memory");
    free(pdata->ghostvals);
    free(pdata->vals);
    free(pdata);
    return NULL;
  }

  for(int i=0; i < NEIGHBOR_TYPE_END; ++i){
    switch (i/2){
    case 0:
      pdata->ghostvals[i] = malloc(sizeof(double) * pnumnodesy*pnumnodesz);
      break;
    case 1:
      pdata->ghostvals[i] = malloc(sizeof(double) * pnumnodesx*pnumnodesz);
    case 2:
      pdata->ghostvals[i] = malloc(sizeof(double) * pnumnodesx*pnumnodesy);
      break;
    default:
      break;
    }
    if(pdata->ghostvals[i] == NULL){
      for(int j = i; j >= 0; --j){
        free(pdata->ghostvals[j]);
      }
      free(pdata->ghostvals);
      free(pdata->vals);
      free(pdata);
      return NULL;
    }
  }

  pdata->grid = *grid;

  return pdata;
}

/******************************************************************************
 * Data file functions                                                        *
 ******************************************************************************/
FILE *create_datafile(grid_t grid, char *filename) {
  if (filename == NULL) {
    DEBUG_PRINT("Invalid NULL filename");
    return NULL;
  }

  FILE *fp;
  if ((fp = fopen(filename, "wb")) == NULL) {
    DEBUG_PRINTF("Failed to open file '%s'", filename);
    return NULL;
  }

  if (fwrite(&grid.numnodesx, sizeof(int), 1, fp) != 1 ||
      fwrite(&grid.numnodesy, sizeof(int), 1, fp) != 1 ||
      fwrite(&grid.numnodesz, sizeof(int), 1, fp) != 1 ||
      fwrite(&grid.xmin, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.xmax, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.ymin, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.ymax, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.zmin, sizeof(double), 1, fp) != 1 ||
      fwrite(&grid.zmax, sizeof(double), 1, fp) != 1) {

    DEBUG_PRINTF("Failed to write header of file '%s'", filename);
    fclose(fp);
    return NULL;
  }

  return fp;
}

FILE *open_datafile(grid_t *grid, int *numsteps, char *filename) {
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
  size_t file_size = ftell(fp);
  rewind(fp);

  if (fread(&grid->numnodesx, sizeof(int), 1, fp) != 1 ||
      fread(&grid->numnodesy, sizeof(int), 1, fp) != 1 ||
      fread(&grid->numnodesz, sizeof(int), 1, fp) != 1 ||
      fread(&grid->xmin, sizeof(double), 1, fp) != 1 ||
      fread(&grid->xmax, sizeof(double), 1, fp) != 1 ||
      fread(&grid->ymin, sizeof(double), 1, fp) != 1 ||
      fread(&grid->ymax, sizeof(double), 1, fp) != 1 ||
      fread(&grid->zmin, sizeof(double), 1, fp) != 1 ||
      fread(&grid->zmax, sizeof(double), 1, fp) != 1) {

    DEBUG_PRINTF("Failed to read header of file '%s'", filename);
    fclose(fp);
    return NULL;
  }

  size_t numnodestot =
      (size_t)grid->numnodesx * grid->numnodesy * grid->numnodesz;

  size_t values_size = numnodestot * sizeof(double);
  size_t stepindex_size = sizeof(int);
  size_t timestamp_size = sizeof(double);
  size_t header_size = 6 * sizeof(double) + 3 * sizeof(int);

  size_t onetimestep_size = values_size + stepindex_size + timestamp_size;
  size_t alltimestep_size = file_size - header_size;

  if (alltimestep_size % onetimestep_size != 0) {
    DEBUG_PRINTF("Data size is inconsistent with number of nodes (%lu, %lu)",
                 alltimestep_size, onetimestep_size);

    fclose(fp);
    return NULL;
  }

  if (numsteps != NULL) {
    *numsteps = (alltimestep_size / onetimestep_size);
  }

  return fp;
}

data_t *read_data(FILE *fp, grid_t *grid, int *step, double *time) {
  if (fp == NULL) {
    DEBUG_PRINT("Invalid NULL file pointer");
    return NULL;
  }

  double ltime;
  int lstep;

  size_t numnodes = NUMNODESTOT(*grid);

  data_t *data;
  if ((data = allocate_data(grid)) == NULL) {
    DEBUG_PRINT("Failed to allocate data");
    return NULL;
  }

  if (fread(&lstep, sizeof(int), 1, fp) != 1 ||
      fread(&ltime, sizeof(double), 1, fp) != 1 ||
      fread(data->vals, sizeof(double), numnodes, fp) != numnodes) {
    DEBUG_PRINT("Failed to read data");
    free(data);
    return NULL;
  }

  if (step != NULL)
    *step = lstep;
  if (time != NULL)
    *time = ltime;

  return data;
}

int write_data(FILE *fp, data_t *data, int step, double time) {
  if (fp == NULL || data == NULL || data->vals == NULL) {
    DEBUG_PRINT("Invalid NULL data or file pointer");
    return 1;
  }

  size_t numnodes = NUMNODESTOT(data->grid);
  if (numnodes <= 0) {
    DEBUG_PRINTF("Invalid number of nodes (%lu)", numnodes);
    return 1;
  }

  if (fwrite(&step, sizeof(int), 1, fp) != 1 ||
      fwrite(&time, sizeof(double), 1, fp) != 1 ||
      fwrite(data->vals, sizeof(double), numnodes, fp) != numnodes) {
    DEBUG_PRINT("Failed to write data");
    return 1;
  }

  return 0;
}

/******************************************************************************
 * Output file functions                                                      *
 ******************************************************************************/
//TODO: Needs to be modified so only one process writes the file
int write_output(output_t *output, data_t *data, int step, double time) {
  if (output == NULL || data == NULL) {
    DEBUG_PRINT("NULL pointer passed as argument");
    return 1;
  }

  output_type_t type = output->type;

  if (type == ALL) {
    return write_data(output->fp, data, step, time);
  }

  int m, n, p;
  closest_index(&data->grid, output->posx, output->posy, output->posz, &m, &n,
                &p);

  int startm = (type == CUTX || type == POINT) ? m : 0;
  int startn = (type == CUTY || type == POINT) ? n : 0;
  int startp = (type == CUTZ || type == POINT) ? p : 0;

  int endm = (type == CUTX || type == POINT) ? m + 1 : NUMNODESX(data);
  int endn = (type == CUTY || type == POINT) ? n + 1 : NUMNODESY(data);
  int endp = (type == CUTZ || type == POINT) ? p + 1 : NUMNODESZ(data);

  data_t *tmpdata = allocate_data(&output->grid);

  for (p = startp; p < endp; p++) {
    for (n = startn; n < endn; n++) {
      for (m = startm; m < endm; m++) {
        int tmpm = m - startm;
        int tmpn = n - startn;
        int tmpp = p - startp;

        SETVALUE(tmpdata, tmpm, tmpn, tmpp, GETVALUE(data, m, n, p));
      }
    }
  }

  int writeok = (write_data(output->fp, tmpdata, step, time) == 0);

  free(tmpdata->vals);
  free(tmpdata);

  if (writeok == 0) {
    DEBUG_PRINT("Failed to write output data");
    return 1;
  }

  return 0;
}

//TODO: Needs to be modified so only one process opens the file
int open_outputfile(output_t *output, grid_t *simgrid) {
  if (output == NULL || simgrid == NULL) {
    DEBUG_PRINT("Invalid NULL pointer in argment");
    return 1;
  }

  grid_t grid;

  output_type_t type = output->type;

  grid.numnodesx = (type == POINT || type == CUTX) ? 1 : simgrid->numnodesx;
  grid.numnodesy = (type == POINT || type == CUTY) ? 1 : simgrid->numnodesy;
  grid.numnodesz = (type == POINT || type == CUTZ) ? 1 : simgrid->numnodesz;

  grid.xmin = (type == POINT || type == CUTX) ? output->posx : simgrid->xmin;
  grid.xmax = (type == POINT || type == CUTX) ? output->posx : simgrid->xmax;

  grid.ymin = (type == POINT || type == CUTY) ? output->posy : simgrid->ymin;
  grid.ymax = (type == POINT || type == CUTY) ? output->posy : simgrid->ymax;

  grid.zmin = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmin;
  grid.zmax = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmax;

  FILE *fp;
  if ((fp = create_datafile(grid, output->filename)) == NULL) {
    DEBUG_PRINTF("Failed to open output file: '%s'", output->filename);
    return 1;
  }

  output->grid = grid;
  output->fp = fp;

  return 0;
}

/******************************************************************************
 * Parameter file functions                                                   *
 ******************************************************************************/
int read_audiosource(char *filename, source_t *source) {
  FILE *fp;
  if ((fp = fopen(filename, "rb")) == NULL) {
    DEBUG_PRINTF("Could not open source file '%s'", filename);
    return 1;
  }

  fseek(fp, 0, SEEK_END);
  size_t filesize = ftell(fp);
  rewind(fp);

  int numsamples = (filesize - sizeof(int)) / sizeof(double);

  int sampling;
  if (fread(&sampling, sizeof(int), 1, fp) != 1) {
    DEBUG_PRINT("Failed to read source data");
    fclose(fp);
    return 1;
  }

  double *data;
  if ((data = malloc(numsamples * sizeof(double))) == NULL) {
    DEBUG_PRINT("Failed to allocate memory for source data");
    return 1;
  }

  int readok = (fread(data, sizeof(double), numsamples, fp) == numsamples);

  fclose(fp);

  if (readok == 0) {
    DEBUG_PRINT("Failed to read source data");
    return 1;
  }

  source->data = data;
  source->numsamples = numsamples;
  source->sampling = sampling;

  return 0;
}

int read_outputparam(FILE *fp, output_t *output) {
  if (fp == NULL || output == NULL) {
    DEBUG_PRINT("NULL passed as argement");
    return 1;
  }

  char typekeyword[BUFSZ_SMALL];
  char sourcekeyword[BUFSZ_SMALL];
  char filename[BUFSZ_LARGE];

  double posxyz[3] = {0.0, 0.0, 0.0};

  if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1 ||
      fscanf(fp, BUFFMT_SMALL, sourcekeyword) != 1 ||
      fscanf(fp, BUFFMT_LARGE, filename) != 1) {

    DEBUG_PRINT("Failed to read an output parameter");
    return 1;
  }

  output_type_t type = CUTX;
  while (type < OUTPUT_TYPE_END &&
         strcmp(output_type_keywords[type], typekeyword) != 0) {
    type++;
  }

  if (type == OUTPUT_TYPE_END) {
    DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
    return 1;
  }

  output_source_t source = PRESSURE;
  while (source < OUTPUT_SOURCE_END &&
         strcmp(output_source_keywords[source], sourcekeyword) != 0) {
    source++;
  }

  if (source == OUTPUT_SOURCE_END) {
    DEBUG_PRINTF("Invalid keyword: '%s'", sourcekeyword);
    return 1;
  }

  int readok = 1;
  switch (type) {
  case CUTX:
    readok = (fscanf(fp, "%lf", &posxyz[0]) == 1);
    break;
  case CUTY:
    readok = (fscanf(fp, "%lf", &posxyz[1]) == 1);
    break;
  case CUTZ:
    readok = (fscanf(fp, "%lf", &posxyz[2]) == 1);
    break;
  case ALL:
    break;

  case POINT:
    readok =
        (fscanf(fp, "%lf %lf %lf", &posxyz[0], &posxyz[1], &posxyz[2]) == 3);
    break;

  default:
    break;
  }

  if (readok == 0) {
    DEBUG_PRINT("Failed to read an output parameter");
    return 1;
  }

  output->filename = copy_string(filename);
  output->type = type;
  output->source = source;
  output->posx = posxyz[0];
  output->posy = posxyz[1];
  output->posz = posxyz[2];

  return 0;
}

int read_sourceparam(FILE *fp, source_t *source) {
  char typekeyword[BUFSZ_SMALL];
  char filename[BUFSZ_LARGE];

  double freq, posx, posy, posz;

  if (fscanf(fp, BUFFMT_SMALL, typekeyword) != 1) {
    DEBUG_PRINT("Failed to read the source parameter");
    return 1;
  }

  source_type_t type = SINE;
  while (type < SOURCE_TYPE_END &&
         strcmp(source_type_keywords[type], typekeyword) != 0) {
    type++;
  }

  if (type == SOURCE_TYPE_END) {
    DEBUG_PRINTF("Invalid keyword: '%s'", typekeyword);
    return 1;
  }

  int readok = 1;
  switch (type) {
  case SINE:
    readok = (fscanf(fp, "%lf", &freq) == 1);
    break;
  case AUDIO:
    readok = (fscanf(fp, BUFFMT_LARGE, filename) == 1);
    break;

  default:
    break;
  }

  if (readok == 0 || fscanf(fp, "%lf %lf %lf", &posx, &posy, &posz) != 3) {
    DEBUG_PRINT("Failed to read the source parameter");
    return 1;
  }

  switch (type) {
  case AUDIO:
    read_audiosource(filename, source);
    break;
  case SINE: {
    if ((source->data = malloc(sizeof(double))) == NULL) {
      DEBUG_PRINT("Failed to allocate memory");
      return 1;
    }

    source->data[0] = freq;
    source->numsamples = 1;

    break;
  }

  default:
    break;
  }

  source->type = type;
  source->posx = posx;
  source->posy = posy;
  source->posz = posz;

  return 0;
}

int read_paramfile(parameters_t *params, const char *filename) {
  if (params == NULL || filename == NULL) {
    DEBUG_PRINT("Invalid print_out params or filename");
    return 1;
  }

  int outrate, numoutputs = 0;

  double dx, dt, maxt;

  char cin_filename[BUFSZ_LARGE];
  char rhoin_filename[BUFSZ_LARGE];

  source_t source;
  output_t *outputs = NULL;

  if ((outputs = malloc(sizeof(output_t) * MAX_OUTPUTS)) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  FILE *fp;
  if ((fp = fopen(filename, "r")) == NULL) {
    DEBUG_PRINTF("Could not open parameter file '%s'", filename);
    return 1;
  }

  int readok =
      ((fscanf(fp, "%lf", &dx) == 1) && (fscanf(fp, "%lf", &dt) == 1) &&
       (fscanf(fp, "%lf", &maxt) == 1) && (fscanf(fp, "%d", &outrate) == 1) &&
       (fscanf(fp, BUFFMT_LARGE, cin_filename) == 1) &&
       (fscanf(fp, BUFFMT_LARGE, rhoin_filename) == 1));

  readok = (readok != 0 && read_sourceparam(fp, &source) == 0 &&fscanf(fp, " ") == 0);

  while (readok != 0 && numoutputs < MAX_OUTPUTS && feof(fp) == 0) {
    readok = (read_outputparam(fp, &outputs[numoutputs++]) == 0 && fscanf(fp, " ") == 0);
  }

  fclose(fp);

  if (readok == 0) {
    DEBUG_PRINT("Failed to read parameter file");
    free(outputs);
    return 1;
  }
  if (numoutputs == 0) {
    free(outputs);
    outputs = NULL;

  }else if((outputs = realloc(outputs, sizeof(output_t) * numoutputs)) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  params->dx = dx;
  params->dt = dt;
  params->maxt = maxt;
  params->outrate = outrate;
  params->cin_filename = copy_string(cin_filename);
  params->rhoin_filename = copy_string(rhoin_filename);
  params->source = source;
  params->numoutputs = numoutputs;
  params->outputs = outputs;

  return 0;
}

/******************************************************************************
 * Simulation-related functions                                               *
 ******************************************************************************/
void apply_source(process_simulation_data_t *simdata, int step) {
  source_t *source = &simdata->params.source;

  double posx = source->posx;
  double posy = source->posy;
  double posz = source->posz;

  double t = step * simdata->params.dt;

  int m, n, p;
  closest_index(&simdata->global_grid, posx, posy, posz, &m, &n, &p);

  if(m < STARTM(simdata->pold) || m > ENDM(simdata->pold) ||
    n < STARTN(simdata->pold) || n > ENDN(simdata->pold) ||
    p < STARTP(simdata->pold) || p > ENDP(simdata->pold)){
    return;
  }

  if (source->type == SINE) {
    double freq = source->data[0];

    process_setvalue(simdata->pold, m, n, p, sin(2 * M_PI * freq * t));

  } else if (source->type == AUDIO) {
    int sample = MIN((int)(t * source->sampling), source->numsamples);

    process_setvalue(simdata->pold, m, n, p, simdata->params.source.data[sample]);
  }
}

int interpolate_inputmaps(process_simulation_data_t *psimdata, process_grid_t *psimgrid,
                          data_t *cin, data_t *rhoin) {
  if (psimdata == NULL || cin == NULL || rhoin == NULL) {
    DEBUG_PRINT("Invalid NULL simdata or cin or rhoin");
    return 1;
  }

  if ((psimdata->c = allocate_pdata(psimgrid)) == NULL ||
      (psimdata->rho = allocate_pdata(psimgrid)) == NULL ||
      (psimdata->rhohalf = allocate_pdata(psimgrid)) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  double dx = psimdata->params.dx;
  double dxd2 = psimdata->params.dx / 2;

  for (int p = psimgrid->startp; p <= psimgrid->endp; p++) {
    for (int n = psimgrid->startn; n <= psimgrid->endn; n++) {
      for (int m = psimgrid->startm; m <= psimgrid->endm; m++) {

        double x = m * dx;
        double y = n * dx;
        double z = p * dx;

        process_setvalue(psimdata->c, m, n, p, trilinear_interpolation(cin, x,y,z));
        process_setvalue(psimdata->rho, m, n, p, trilinear_interpolation(rhoin, x,y,z));

        x += dxd2;
        y += dxd2;
        z += dxd2;

        process_setvalue(psimdata->rhohalf, m, n, p, trilinear_interpolation(rhoin, x,y,z));
      }
    }
  }

  return 0;
}

//TODO: Needs to accept a process_simulation_data_t*
//TODO: processes need to send their old velocities to their RIGHT,BACK,UP neighbor
void update_pressure(simulation_data_t *simdata) {
  const double dtdx = simdata->params.dt / simdata->params.dx;

  const int numnodesx = NUMNODESX(simdata->pold);
  const int numnodesy = NUMNODESY(simdata->pold);
  const int numnodesz = NUMNODESZ(simdata->pold);

  for (int p = 0; p < numnodesz; p++) {
    for (int n = 0; n < numnodesy; n++) {
      for (int m = 0; m < numnodesx; m++) {
        double c = GETVALUE(simdata->c, m, n, p);
        double rho = GETVALUE(simdata->rho, m, n, p);

        double rhoc2dtdx = rho * c * c * dtdx;

        double dvx = GETVALUE(simdata->vxold, m, n, p);
        double dvy = GETVALUE(simdata->vyold, m, n, p);
        double dvz = GETVALUE(simdata->vzold, m, n, p);

        dvx -= m > 0 ? GETVALUE(simdata->vxold, m - 1, n, p) : 0.0;
        dvy -= n > 0 ? GETVALUE(simdata->vyold, m, n - 1, p) : 0.0;
        dvz -= p > 0 ? GETVALUE(simdata->vzold, m, n, p - 1) : 0.0;

        double prev_p = GETVALUE(simdata->pold, m, n, p);

        SETVALUE(simdata->pnew, m, n, p,
                 prev_p - rhoc2dtdx * (dvx + dvy + dvz));
      }
    }
  }
}

//TODO: Needs to accept a process_simulation_data_t*
//TODO: processes need to send their new pressures to theur LEFT,FRONT,DOWN neighbor
void update_velocities(simulation_data_t *simdata) {
  const double dtdx = simdata->params.dt / simdata->params.dx;

  const int numnodesx = NUMNODESX(simdata->vxold);
  const int numnodesy = NUMNODESY(simdata->vxold);
  const int numnodesz = NUMNODESZ(simdata->vxold);

  for (int p = 0; p < numnodesz; p++) {
    for (int n = 0; n < numnodesy; n++) {
      for (int m = 0; m < numnodesx; m++) {
        int mp1 = MIN(numnodesx - 1, m + 1);
        int np1 = MIN(numnodesy - 1, n + 1);
        int pp1 = MIN(numnodesz - 1, p + 1);

        double dtdxrho = dtdx / GETVALUE(simdata->rhohalf, m, n, p);

        double p_mnq = GETVALUE(simdata->pnew, m, n, p);

        double dpx = GETVALUE(simdata->pnew, mp1, n, p) - p_mnq;
        double dpy = GETVALUE(simdata->pnew, m, np1, p) - p_mnq;
        double dpz = GETVALUE(simdata->pnew, m, n, pp1) - p_mnq;

        double prev_vx = GETVALUE(simdata->vxold, m, n, p);
        double prev_vy = GETVALUE(simdata->vyold, m, n, p);
        double prev_vz = GETVALUE(simdata->vzold, m, n, p);

        SETVALUE(simdata->vxnew, m, n, p, prev_vx - dtdxrho * dpx);
        SETVALUE(simdata->vynew, m, n, p, prev_vy - dtdxrho * dpy);
        SETVALUE(simdata->vznew, m, n, p, prev_vz - dtdxrho * dpz);
      }
    }
  }
}

void init_simulation(process_simulation_data_t *psimdata, const char *params_filename) {
  if (read_paramfile(&psimdata->params, params_filename) != 0) {
    printf("Failed to read parameters. Aborting...\n\n");
    exit(1);
  }

  grid_t rhoin_grid;
  grid_t cin_grid;
  process_grid_t psim_grid;

  int rho_numstep;
  int c_numstep;

  FILE *rhofp = open_datafile(&rhoin_grid, &rho_numstep, psimdata->params.rhoin_filename);
  FILE *cfp = open_datafile(&cin_grid, &c_numstep, psimdata->params.cin_filename);

  if (rhofp == NULL || rho_numstep <= 0) {
    printf("Failed to open the density map file. Aborting...\n\n");
    exit(1);
  }

  if (cfp == NULL || c_numstep <= 0) {
    printf("Failed to open the speed map file. Aborting...\n\n");
    exit(1);
  }

  if (rhoin_grid.xmin != cin_grid.xmin || rhoin_grid.ymin != cin_grid.ymin ||
      rhoin_grid.zmin != cin_grid.zmin || rhoin_grid.xmax != cin_grid.xmax ||
      rhoin_grid.ymax != cin_grid.ymax || rhoin_grid.zmax != cin_grid.zmax) {
    printf("Grids for the density and speed are not the same. Aborting...\n\n");
    exit(1);
  }

  data_t *rho_map = read_data(rhofp, &rhoin_grid, NULL, NULL);
  data_t *c_map = read_data(cfp, &cin_grid, NULL, NULL);

  if (rho_map == NULL || c_map == NULL) {
    printf("Failed to read data from input maps. Aborting...\n\n");
    exit(1);
  }

  fclose(rhofp);
  fclose(cfp);

  double xmin = rhoin_grid.xmin, xmax = rhoin_grid.xmax, 
    ymin = rhoin_grid.ymin, ymax = rhoin_grid.ymax, 
    zmin = rhoin_grid.zmin, zmax = rhoin_grid.zmax;

  int numnodesx = MAX(floor((xmax - xmin) / psimdata->params.dx), 1),
    numnodesy = MAX(floor((ymax - ymin) / psimdata->params.dx), 1),
    numnodesz = MAX(floor((zmax - zmin) / psimdata->params.dx), 1);

  // Copy global grid into process simulation data
  psimdata->global_grid.xmin = xmin; psimdata->global_grid.xmax = xmax;
  psimdata->global_grid.ymin = ymin; psimdata->global_grid.ymax = ymax;
  psimdata->global_grid.zmin = zmin; psimdata->global_grid.zmax = zmax;

  psimdata->global_grid.numnodesx = numnodesx;
  psimdata->global_grid.numnodesy = numnodesy;
  psimdata->global_grid.numnodesz = numnodesz;

  //Set correct value ranges
  psim_grid.gnumx = numnodesx; psim_grid.gnumy = numnodesy; psim_grid.gnumz = numnodesz;
  psim_grid.startm = numnodesx*coords[0]/dims[0]; psim_grid.endm = numnodesx*(coords[0]+1)/dims[0] - 1;
  psim_grid.startn = numnodesy*coords[1]/dims[1]; psim_grid.endn = numnodesy*(coords[1]+1)/dims[1] - 1;
  psim_grid.startp = numnodesz*coords[2]/dims[2]; psim_grid.endp = numnodesz*(coords[2]+1)/dims[2] - 1;

  if (interpolate_inputmaps(psimdata, &psim_grid, c_map, rho_map) != 0) {
    printf(
        "Error while converting input map to simulation grid. Aborting...\n\n");
    exit(1);
  }

  printf("Rank %4d has subdomain (%3d, %3d) (%3d, %3d) (%3d, %3d) of global grid (%3d, %3d, %3d)\n\tcvalue %10.5lf at (%3d, %3d, %3d)\n",
    cart_rank, psim_grid.startm, psim_grid.endm, psim_grid.startn, psim_grid.endn, psim_grid.startp, psim_grid.endp,
    psim_grid.gnumx, psim_grid.gnumy, psim_grid.gnumz, psimdata->c->vals[0], psim_grid.startm, psim_grid.startn, psim_grid.startp);

#if 0
  if (simdata->params.outrate > 0 && simdata->params.outputs != NULL) {
    for (int i = 0; i < simdata->params.numoutputs; i++) {
      char *outfilei = simdata->params.outputs[i].filename;

      for (int j = 0; j < i; j++) {
        char *outfilej = simdata->params.outputs[j].filename;

        if (strcmp(outfilei, outfilej) == 0) {
          printf("Duplicate output file: '%s'. Aborting...\n\n", outfilei);
          exit(1);
        }
      }
    }

    for (int i = 0; i < simdata->params.numoutputs; i++) {
      output_t *output = &simdata->params.outputs[i];

      if (open_outputfile(output, &sim_grid) != 0) {
        printf("Failed to open output file: '%s'. Aborting...\n\n",
               output->filename);
        exit(1);
      }
    }
  }

  if ((simdata->pold = allocate_data(&sim_grid)) == NULL ||
      (simdata->pnew = allocate_data(&sim_grid)) == NULL ||
      (simdata->vxold = allocate_data(&sim_grid)) == NULL ||
      (simdata->vxnew = allocate_data(&sim_grid)) == NULL ||
      (simdata->vyold = allocate_data(&sim_grid)) == NULL ||
      (simdata->vynew = allocate_data(&sim_grid)) == NULL ||
      (simdata->vzold = allocate_data(&sim_grid)) == NULL ||
      (simdata->vznew = allocate_data(&sim_grid)) == NULL) {
    printf("Failed to allocate memory. Aborting...\n\n");
    exit(1);
  }

  fill_data(simdata->pold, 0.0);
  fill_data(simdata->pnew, 0.0);

  fill_data(simdata->vynew, 0.0);
  fill_data(simdata->vxold, 0.0);
  fill_data(simdata->vynew, 0.0);
  fill_data(simdata->vyold, 0.0);
  fill_data(simdata->vznew, 0.0);
  fill_data(simdata->vzold, 0.0);

  printf("\n");
  printf(" Grid spacing: %g\n", simdata->params.dx);
  printf("  Grid size X: %d\n", sim_grid.numnodesx);
  printf("  Grid size Y: %d\n", sim_grid.numnodesy);
  printf("  Grid size Z: %d\n", sim_grid.numnodesz);
  printf("    Time step: %g\n", simdata->params.dt);
  printf(" Maximum time: %g\n\n", simdata->params.maxt);

  if (simdata->params.outrate > 0 && simdata->params.outputs) {
    int outsampling =
        (int)(1.0 / (simdata->params.outrate * simdata->params.dt));

    printf("     Output rate: every %d step(s)\n", simdata->params.outrate);
    printf(" Output sampling: %d Hz\n\n", outsampling);
    printf(" Output files:\n\n");

    for (int i = 0; i < simdata->params.numoutputs; i++) {
      print_output(&simdata->params.outputs[i]);
    }

    printf("\n");

  } else if (simdata->params.outrate < 0) {
    printf("  Output is disabled (output rate set to 0)\n\n");

  } else {
    printf("  Output is disabled (not output specified)\n\n");
  }

  print_source(&simdata->params.source);

  fflush(stdout);

#endif

  free(rho_map->vals);
  free(rho_map);
  free(c_map->vals);
  free(c_map);
}

//TODO: Needs to be modified so only one process closes the files
void finalize_simulation(process_simulation_data_t *simdata) {
  if (simdata->params.outputs != NULL) {
    for (int i = 0; i < simdata->params.numoutputs; i++) {
      free(simdata->params.outputs[i].filename);

      /*if (simdata->params.outrate > 0) {
        fclose(simdata->params.outputs[i].fp);
      }*/
    }

    free(simdata->params.outputs);
  }

  free(simdata->params.source.data);
  free(simdata->params.cin_filename);
  free(simdata->params.rhoin_filename);

  return; // return early, not everything is malloc yet

  free(simdata->rho->vals);
  free(simdata->rho);
  free(simdata->rhohalf->vals);
  free(simdata->rhohalf);
  free(simdata->c->vals);
  free(simdata->c);

#if 0
  free(simdata->pold->vals);
  free(simdata->pold);
  free(simdata->pnew->vals);
  free(simdata->pnew);

  free(simdata->vxold->vals);
  free(simdata->vxold);
  free(simdata->vxnew->vals);
  free(simdata->vxnew);
  free(simdata->vyold->vals);
  free(simdata->vyold);
  free(simdata->vynew->vals);
  free(simdata->vynew);
  free(simdata->vzold->vals);
  free(simdata->vzold);
  free(simdata->vznew->vals);
  free(simdata->vznew);
#endif
}

void swap_timesteps(process_simulation_data_t *psimdata) {
  process_data_t *tmpp = psimdata->pold;
  process_data_t *tmpvx = psimdata->vxold;
  process_data_t *tmpvy = psimdata->vyold;
  process_data_t *tmpvz = psimdata->vzold;

  psimdata->pold = psimdata->pnew;
  psimdata->pnew = tmpp;
  psimdata->vxold = psimdata->vxnew;
  psimdata->vxnew = tmpvx;
  psimdata->vyold = psimdata->vynew;
  psimdata->vynew = tmpvy;
  psimdata->vzold = psimdata->vznew;
  psimdata->vznew = tmpvz;
}