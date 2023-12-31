#include <mpi.h>
#include <stdio.h>

#include "fdtd.h"

// Global topology values that once set never change
MPI_Comm cart_comm;
int cart_rank;
int dims[3];
int coords[3];
int neighbors[6];

#pragma omp declare mapper(source_t source) \
              map(source) \
              map(source.data[0:source.numsamples])

#pragma omp declare mapper(process_data_t data) \
              map(data) \
              map(data.vals[0:PNUMNODESTOT(&data)])

#pragma omp declare mapper(process_simulation_data_t simdata) \
              map(simdata) \
              map(simdata.c[0:1]) map(simdata.rho[0:1]) map(simdata.rhohalf[0:1]) \
              map(simdata.pold[0:1]) map(simdata.pold->ghostvals[0:NEIGHBOR_TYPE_END]) \
              map(simdata.pold->ghostvals[RIGHT][0:PNUMNODESY(simdata.pold) * PNUMNODESZ(simdata.pold)]) \
              map(simdata.pold->ghostvals[BACK ][0:PNUMNODESX(simdata.pold) * PNUMNODESZ(simdata.pold)]) \
              map(simdata.pold->ghostvals[UP   ][0:PNUMNODESX(simdata.pold) * PNUMNODESY(simdata.pold)]) \
              map(simdata.pnew[0:1]) map(simdata.pnew->ghostvals[0:NEIGHBOR_TYPE_END]) \
              map(simdata.pnew->ghostvals[RIGHT][0:PNUMNODESY(simdata.pnew) * PNUMNODESZ(simdata.pnew)]) \
              map(simdata.pnew->ghostvals[BACK ][0:PNUMNODESX(simdata.pnew) * PNUMNODESZ(simdata.pnew)]) \
              map(simdata.pnew->ghostvals[UP   ][0:PNUMNODESX(simdata.pnew) * PNUMNODESY(simdata.pnew)]) \
              map(simdata.vxold[0:1]) map(simdata.vxold->ghostvals[0:NEIGHBOR_TYPE_END])\
              map(simdata.vxold->ghostvals[LEFT ][0:PNUMNODESY(simdata.vxold) * PNUMNODESZ(simdata.vxold)]) \
              map(simdata.vxnew[0:1]) map(simdata.vxnew->ghostvals[0:NEIGHBOR_TYPE_END])\
              map(simdata.vxnew->ghostvals[LEFT ][0:PNUMNODESY(simdata.vxnew) * PNUMNODESZ(simdata.vxnew)]) \
              map(simdata.vyold[0:1]) map(simdata.vyold->ghostvals[0:NEIGHBOR_TYPE_END])\
              map(simdata.vyold->ghostvals[FRONT][0:PNUMNODESX(simdata.vyold) * PNUMNODESZ(simdata.vyold)]) \
              map(simdata.vynew[0:1]) map(simdata.vynew->ghostvals[0:NEIGHBOR_TYPE_END])\
              map(simdata.vynew->ghostvals[FRONT][0:PNUMNODESX(simdata.vynew) * PNUMNODESZ(simdata.vynew)]) \
              map(simdata.vzold[0:1]) map(simdata.vzold->ghostvals[0:NEIGHBOR_TYPE_END])\
              map(simdata.vzold->ghostvals[DOWN ][0:PNUMNODESX(simdata.vzold) * PNUMNODESY(simdata.vzold)]) \
              map(simdata.vznew[0:1]) map(simdata.vznew->ghostvals[0:NEIGHBOR_TYPE_END])\
              map(simdata.vznew->ghostvals[DOWN ][0:PNUMNODESX(simdata.vznew) * PNUMNODESY(simdata.vznew)]) \
              map(simdata.buffer_vx[0:PNUMNODESY(simdata.pold) * PNUMNODESZ(simdata.pold)]) \
              map(simdata.buffer_px[0:PNUMNODESY(simdata.pold) * PNUMNODESZ(simdata.pold)]) \
              map(simdata.buffer_vy[0:PNUMNODESX(simdata.pold) * PNUMNODESZ(simdata.pold)]) \
              map(simdata.buffer_py[0:PNUMNODESX(simdata.pold) * PNUMNODESZ(simdata.pold)]) \
              map(simdata.buffer_vz[0:PNUMNODESX(simdata.pold) * PNUMNODESY(simdata.pold)]) \
              map(simdata.buffer_pz[0:PNUMNODESX(simdata.pold) * PNUMNODESY(simdata.pold)])

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
  MPI_Comm_size(cart_comm, &world_size);
  MPI_Comm_rank(cart_comm, &cart_rank);

  MPI_Cart_coords(cart_comm, cart_rank, 3, coords);
  MPI_Cart_shift(cart_comm, 0, 1, 
    &neighbors[LEFT], &neighbors[RIGHT]);
  MPI_Cart_shift(cart_comm, 1, 1, 
    &neighbors[FRONT], &neighbors[BACK]);
  MPI_Cart_shift(cart_comm, 2, 1, 
		&neighbors[DOWN], &neighbors[UP]);


  process_simulation_data_t psimdata;
  init_simulation(&psimdata, argv[1]);

  int numtimesteps = floor(psimdata.params.maxt / psimdata.params.dt);

  double start = GET_TIME();
  double apply_time, output_time, updateP_time, updateV_time, swap_time;
  double t1, t2;
  #pragma omp target data map(to:psimdata)
  for (int tstep = 0; tstep <= numtimesteps; tstep++) {
    t1 = GET_TIME(); 
    apply_source(&psimdata, tstep);
    t2 = GET_TIME();
    apply_time += t2 - t1;
    // Need to collect data from all the different processes
    if (psimdata.params.outrate > 0 && (tstep % psimdata.params.outrate) == 0) {
      int send_size = PNUMNODESX(psimdata.pold)*PNUMNODESY(psimdata.pold)*PNUMNODESZ(psimdata.pold);
      for (int i = 0; i < psimdata.params.numoutputs; i++) {
        process_data_t *output_data = NULL;

        switch (psimdata.params.outputs[i].source) {
          case PRESSURE:
            #pragma omp target update from(psimdata.pold[0:1])
            output_data = psimdata.pold;
            break;
          case VELOCITYX:
            #pragma omp target update from(psimdata.vxold[0:1])
            output_data = psimdata.vxold;
            break;
          case VELOCITYY:
            #pragma omp target update from(psimdata.vyold[0:1])
            output_data = psimdata.vyold;
            break;
          case VELOCITYZ:
            #pragma omp target update from(psimdata.vzold[0:1])
            output_data = psimdata.vzold;
            break;

          default:
            break;
        }

        double time = tstep * psimdata.params.dt;
        write_output(&psimdata.params.outputs[i], output_data, &psimdata.global_grid ,tstep, time);
      }
    }
    t1 = GET_TIME();
    output_time += t1 - t2;

    if(cart_rank == 0){
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
    }
    t1 = GET_TIME();
    update_pressure(&psimdata);
    t2 = GET_TIME();
    updateP_time += t2 - t1;
    update_velocities(&psimdata);
    t1 = GET_TIME();
    updateV_time += t1 - t2;
    swap_timesteps(&psimdata);
    t2 = GET_TIME();
    swap_time += t2 - t1;
  }

  printf("%d : %10.3lf | %10.3lf | %10.3lf | %10.3lf | %10.3lf\n", 
    cart_rank, apply_time, output_time, updateP_time, updateV_time, swap_time);

  if(cart_rank == 0){
    double elapsed = GET_TIME() - start;
    double numupdates =
        (double)NUMNODESTOT(psimdata.global_grid) * (numtimesteps + 1);
    double updatespers = numupdates / elapsed / 1e6;
    printf("\nElapsed %.6lf seconds (%.3lf Mupdates/s)\n\n", elapsed, updatespers);
    fflush(stdout);
  }

  finalize_simulation(&psimdata);

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

inline double process_getvalue(process_data_t *pdata, int m, int n, int p) {
  
  int mbar = m - STARTM(pdata), nbar = n - STARTN(pdata), pbar = p - STARTP(pdata);
  int pnumnodesx = PNUMNODESX(pdata), pnumnodesy = PNUMNODESY(pdata), pnumnodesz = PNUMNODESZ(pdata);

  if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
     nbar >= 0 && nbar <= pnumnodesy - 1 && 
     pbar >= 0 && pbar <= pnumnodesz - 1){
    return pdata->vals[pbar*pnumnodesy*pnumnodesx + nbar*pnumnodesx + mbar];
  }
  else if(nbar >= 0 && nbar <= pnumnodesy - 1 &&
          pbar >= 0 && pbar <= pnumnodesz - 1) {
    return pdata->ghostvals[(mbar == -1)?LEFT:RIGHT][pbar*pnumnodesy + nbar];
  }
  else if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
          pbar >= 0 && pbar <= pnumnodesz - 1) {
    return pdata->ghostvals[(nbar == -1)?FRONT:BACK][pbar*pnumnodesx + mbar];
  }
  else if(mbar >= 0 && mbar <= pnumnodesx - 1 &&
          nbar >= 0 && nbar <= pnumnodesy - 1) {
    return pdata->ghostvals[(pbar == -1)?DOWN:UP][nbar*pnumnodesx + mbar];
  }
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

void fill_data(process_data_t *pdata, double value) {
  if (pdata == NULL) {
    DEBUG_PRINT("Invalid NULL data");
    return;
  }

  #pragma omp parallel for collapse(2)
  for (int pbar = 0; pbar < PNUMNODESZ(pdata); pbar++) {
    for (int nbar = 0; nbar < PNUMNODESY(pdata); nbar++) {
      for (int mbar = 0; mbar < PNUMNODESX(pdata); mbar++) {
        PROCESS_SETVALUE_INSIDE(pdata, mbar, nbar, pbar, value);
      }
    }
  }
}

process_data_t *allocate_pdata(process_grid_t *grid, int malloc_ghost_flags) {
  int pnumnodesx = grid->lm.n, 
      pnumnodesy = grid->ln.n, 
      pnumnodesz = grid->lp.n;
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

  pdata->malloc_ghost_flags = malloc_ghost_flags;

  for(int i=0; i < NEIGHBOR_TYPE_END; ++i){
    // Only allocate asked ghost cells
    if(!((1 << i) & malloc_ghost_flags)) continue;
    switch (i/2){
    case 0:
      pdata->ghostvals[i] = malloc(sizeof(double) * pnumnodesy*pnumnodesz);
      break;
    case 1:
      pdata->ghostvals[i] = malloc(sizeof(double) * pnumnodesx*pnumnodesz);
      break;
    case 2:
      pdata->ghostvals[i] = malloc(sizeof(double) * pnumnodesx*pnumnodesy);
      break;
    default:
      break;
    }
    if(pdata->ghostvals[i] == NULL){
      for(int j = i; j >= 0; --j){
        if (malloc_ghost_flags & (1 << j)) free(pdata->ghostvals[j]);
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

void free_pdata(process_data_t *pdata){
  if(pdata == NULL) return;

  for(int i = 0; i < NEIGHBOR_TYPE_END; ++i){
    if (pdata->malloc_ghost_flags & (1 << i)) {
      free(pdata->ghostvals[i]);
    }
  }
  free(pdata->ghostvals);
  free(pdata->vals);
  free(pdata);
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
int write_output(output_t *output, process_data_t *data, grid_t* global_grid, int step, double time) {
  if (output == NULL || data == NULL) {
    DEBUG_PRINT("NULL pointer passed as argument");
    return 1;
  }

  if(output->fp == NULL){
    //This process should not save anything
    return 0;
  }

  output_type_t type = output->type;
  
  /*
  if (type == ALL) {
    return write_data(output->fp, data, step, time);
  }
  */

  int m, n, p;
  closest_index(global_grid, output->posx, output->posy, output->posz, &m, &n,
                &p);

  int mbar = m - data->grid.lm.start,
    nbar = n - data->grid.ln.start,
    pbar = p - data->grid.lp.start;

  int startm = (type == CUTX || type == POINT) ? mbar : 0;
  int startn = (type == CUTY || type == POINT) ? nbar : 0;
  int startp = (type == CUTZ || type == POINT) ? pbar : 0;

  int endm = (type == CUTX || type == POINT) ? mbar + 1 : PNUMNODESX(data);
  int endn = (type == CUTY || type == POINT) ? nbar + 1 : PNUMNODESY(data);
  int endp = (type == CUTZ || type == POINT) ? pbar + 1 : PNUMNODESZ(data);

  #pragma omp target update from(data->vals[PINDEX3D(data, mbar, nbar, pbar):1])

  data_t *tmpdata = allocate_data(&output->grid);

  for (p = startp; p < endp; p++) {
    for (n = startn; n < endn; n++) {
      for (m = startm; m < endm; m++) {
        int tmpm = m - startm;
        int tmpn = n - startn;
        int tmpp = p - startp;

        SETVALUE(tmpdata, tmpm, tmpn, tmpp, PROCESS_GETVALUE_INSIDE(data, m, n, p));
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

int open_outputfile(output_t *output, process_grid_t *simgrid, grid_t* global_grid) {
  if (output == NULL || simgrid == NULL) {
    DEBUG_PRINT("Invalid NULL pointer in argment");
    return 1;
  }

  grid_t grid;

  output_type_t type = output->type;

  grid.numnodesx = (type == POINT || type == CUTX) ? 1 : simgrid->lm.n;
  grid.numnodesy = (type == POINT || type == CUTY) ? 1 : simgrid->ln.n;
  grid.numnodesz = (type == POINT || type == CUTZ) ? 1 : simgrid->lp.n;

  grid.xmin = (type == POINT || type == CUTX) ? output->posx : simgrid->xmin;
  grid.xmax = (type == POINT || type == CUTX) ? output->posx : simgrid->xmax;

  grid.ymin = (type == POINT || type == CUTY) ? output->posy : simgrid->ymin;
  grid.ymax = (type == POINT || type == CUTY) ? output->posy : simgrid->ymax;

  grid.zmin = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmin;
  grid.zmax = (type == POINT || type == CUTZ) ? output->posz : simgrid->zmax;

  int m, n, p;
  closest_index(global_grid, output->posx, output->posy, output->posz, &m, &n, &p);

  char buffer_filename[BUFSZ_LARGE];
  snprintf(buffer_filename, BUFSZ_LARGE, "%d_%d_%d_%s", 
    coords[0], coords[1], coords[2], output->filename);

  FILE *fp = NULL;
  if ((fp = create_datafile(grid, buffer_filename)) == NULL) {
    DEBUG_PRINTF("Failed to open output file: '%s'", buffer_filename);
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

  if (type != ALL){
    DEBUG_PRINT("Limited functionality, only ALL output type accepted");
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
void apply_source(process_simulation_data_t *psimdata, int step) {
  source_t *source = &psimdata->params.source;

  double posx = source->posx;
  double posy = source->posy;
  double posz = source->posz;

  double t = step * psimdata->params.dt;

  int m, n, p;
  closest_index(&psimdata->global_grid, posx, posy, posz, &m, &n, &p);

  register process_data_t *pold = psimdata->pold;
  const register int mbar = m - STARTM(pold), nbar = n - STARTN(pold), pbar = p - STARTP(pold);

  if(mbar < 0 || mbar > PNUMNODESX(pold)-1 ||
     nbar < 0 || nbar > PNUMNODESY(pold)-1 ||
     pbar < 0 || pbar > PNUMNODESZ(pold)-1){
    return;
  }

  if (source->type == SINE) {
    double freq = source->data[0];

    PROCESS_SETVALUE_INSIDE(pold, mbar, nbar, pbar, sin(2 * M_PI * freq * t));

  } else if (source->type == AUDIO) {
    int sample = MIN((int)(t * source->sampling), source->numsamples);

    PROCESS_SETVALUE_INSIDE(pold, mbar, nbar, pbar, psimdata->params.source.data[sample]);
  }

  #pragma omp target update to(psimdata->pold->vals[ PINDEX3D(psimdata->pold, mbar,nbar,pbar) :1])
}

int interpolate_inputmaps(process_simulation_data_t *psimdata, process_grid_t *psimgrid,
                          data_t *cin, data_t *rhoin) {
  if (psimdata == NULL || cin == NULL || rhoin == NULL) {
    DEBUG_PRINT("Invalid NULL simdata or cin or rhoin");
    return 1;
  }

  DEBUG_PRINT("Entered interpolate_inputmaps");

  if ((psimdata->c = allocate_pdata(psimgrid, 0)) == NULL ||
      (psimdata->rho = allocate_pdata(psimgrid, 0)) == NULL ||
      (psimdata->rhohalf = allocate_pdata(psimgrid, 0)) == NULL) {
    DEBUG_PRINT("Failed to allocate memory");
    return 1;
  }

  DEBUG_PRINT("Created inputs");

  double dx = psimdata->params.dx;
  double dxd2 = psimdata->params.dx / 2;

  for (int pbar = 0; pbar < psimgrid->lp.n; pbar++) {
    for (int nbar = 0; nbar < psimgrid->ln.n; nbar++) {
      for (int mbar = 0; mbar < psimgrid->lm.n; mbar++) {

        double x = (mbar + psimgrid->lm.start) * dx;
        double y = (nbar + psimgrid->ln.start) * dx;
        double z = (pbar + psimgrid->lp.start) * dx;

        PROCESS_SETVALUE_INSIDE(psimdata->c      , mbar, nbar, pbar, trilinear_interpolation(cin, x,y,z));
        PROCESS_SETVALUE_INSIDE(psimdata->rho    , mbar, nbar, pbar, trilinear_interpolation(rhoin, x,y,z));
        PROCESS_SETVALUE_INSIDE(psimdata->rhohalf, mbar, nbar, pbar, trilinear_interpolation(rhoin, x+dxd2, y+dxd2, z+dxd2));
      }
    }
  }

  return 0;
}


void update_pressure(process_simulation_data_t *psimdata) {

  process_data_t *pold = psimdata->pold;
  MPI_Request request_recv[3], request_send[3];

  int m, n, p;
  int mbar, nbar, pbar;
  register const int pnumnodesx = PNUMNODESX(pold), startm = STARTM(pold), endm = ENDM(pold),
                     pnumnodesy = PNUMNODESY(pold), startn = STARTN(pold), endn = ENDN(pold),
                     pnumnodesz = PNUMNODESZ(pold), startp = STARTP(pold), endp = ENDP(pold);

  #pragma omp target teams distribute
  for (pbar = 0; pbar < pnumnodesz; pbar++) {
    #pragma omp parallel for
    for (nbar = 0; nbar < pnumnodesy; nbar++){
      psimdata->buffer_vx[pbar * pnumnodesy + nbar] = PROCESS_GETVALUE_INSIDE(psimdata->vxold, pnumnodesx-1, nbar, pbar);
    }
  }
  #pragma omp target teams distribute
  for (pbar = 0; pbar < pnumnodesz; pbar++) {
    #pragma omp parallel for
    for (mbar = 0; mbar < pnumnodesx; mbar++){
      psimdata->buffer_vy[pbar * pnumnodesx + mbar] = PROCESS_GETVALUE_INSIDE(psimdata->vyold, mbar, pnumnodesy-1, pbar);
    }
  }
  #pragma omp target teams distribute
  for (nbar = 0; nbar < pnumnodesy; nbar++) {
    #pragma omp parallel for
    for (mbar = 0; mbar < pnumnodesx; mbar++){
      psimdata->buffer_vz[nbar * pnumnodesx + mbar] = PROCESS_GETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pnumnodesz-1);
    }
  }

  #pragma omp target update from(psimdata->buffer_vx[0:pnumnodesy*pnumnodesz])
  #pragma omp target update from(psimdata->buffer_vy[0:pnumnodesx*pnumnodesz])
  #pragma omp target update from(psimdata->buffer_vz[0:pnumnodesx*pnumnodesy])

  MPI_Isend(psimdata->buffer_vx, pnumnodesy*pnumnodesz, MPI_DOUBLE, neighbors[RIGHT], SEND_X, cart_comm, &request_send[0]);
  MPI_Isend(psimdata->buffer_vy, pnumnodesx*pnumnodesz, MPI_DOUBLE, neighbors[BACK ], SEND_Y, cart_comm, &request_send[1]);
  MPI_Isend(psimdata->buffer_vz, pnumnodesx*pnumnodesy, MPI_DOUBLE, neighbors[UP   ], SEND_Z, cart_comm, &request_send[2]);

  MPI_Irecv(psimdata->vxold->ghostvals[LEFT ], pnumnodesy*pnumnodesz, MPI_DOUBLE, neighbors[LEFT ], SEND_X, cart_comm, &request_recv[0]);
  MPI_Irecv(psimdata->vyold->ghostvals[FRONT], pnumnodesx*pnumnodesz, MPI_DOUBLE, neighbors[FRONT], SEND_Y, cart_comm, &request_recv[1]);
  MPI_Irecv(psimdata->vzold->ghostvals[DOWN ], pnumnodesx*pnumnodesy, MPI_DOUBLE, neighbors[DOWN ], SEND_Z, cart_comm, &request_recv[2]);

  #pragma omp target teams distribute
  for (p = startp + 1; p <= endp; p++) {
    #pragma omp parallel for collapse(2)
    for (n = startn + 1; n <= endn; n++) {
      for (m = startm + 1; m <= endm; m++) {
        update_pressure_routine_inside(psimdata, m, n, p);
      }
    }
  }

  MPI_Waitall(3, request_recv, MPI_STATUSES_IGNORE);

  #pragma omp target update to(psimdata->vxold->ghostvals[LEFT ][0:pnumnodesy*pnumnodesz])
  #pragma omp target update to(psimdata->vyold->ghostvals[FRONT][0:pnumnodesx*pnumnodesz])
  #pragma omp target update to(psimdata->vzold->ghostvals[DOWN ][0:pnumnodesx*pnumnodesy])

  #pragma omp target teams distribute
  for (p = startp; p <= endp; p++) {
    #pragma omp parallel for
    for (n = startn; n <= endn; n++) {
      update_pressure_routine(psimdata, startm, n, p, LEFT);
    }
  }
  #pragma omp target teams distribute
  for (p = startp; p <= endp; p++) {
    #pragma omp parallel for
    for (m = startm + 1; m <= endm; m++) {
      update_pressure_routine(psimdata, m, startn, p, FRONT);
    }
  }
  #pragma omp target teams distribute
  for (n = startn + 1; n <= endn; n++){
    #pragma omp parallel for
    for (m = startm + 1; m <= endm; m++){
      update_pressure_routine(psimdata, m, n, startp, DOWN);
    }
  }

  MPI_Waitall(3, request_send, MPI_STATUSES_IGNORE);
}

void update_pressure_routine_inside(process_simulation_data_t *psimdata, int m, int n, int p) {

  const double dtdx = psimdata->params.dt / psimdata->params.dx;
  const register int mbar = m - STARTM(psimdata->pold), nbar = n - STARTN(psimdata->pold), pbar = p - STARTP(psimdata->pold);

  double c   = PROCESS_GETVALUE_INSIDE(psimdata->c    , mbar, nbar, pbar);
  double rho = PROCESS_GETVALUE_INSIDE(psimdata->rho  , mbar, nbar, pbar);

  double dvx = PROCESS_GETVALUE_INSIDE(psimdata->vxold, mbar, nbar, pbar) - PROCESS_GETVALUE_INSIDE(psimdata->vxold, mbar - 1, nbar, pbar);
  double dvy = PROCESS_GETVALUE_INSIDE(psimdata->vyold, mbar, nbar, pbar) - PROCESS_GETVALUE_INSIDE(psimdata->vyold, mbar, nbar - 1, pbar);
  double dvz = PROCESS_GETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pbar) - PROCESS_GETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pbar - 1);

  double prev_p = PROCESS_GETVALUE_INSIDE(psimdata->pold, mbar, nbar, pbar);

  PROCESS_SETVALUE_INSIDE(psimdata->pnew, mbar, nbar, pbar, prev_p - rho * c * c * dtdx * (dvx + dvy + dvz));
}

void update_pressure_routine(process_simulation_data_t *psimdata, int m, int n, int p, neighbor_t neighbor) {

  const double dtdx = psimdata->params.dt / psimdata->params.dx;
  const register int mbar = m - STARTM(psimdata->pold), nbar = n - STARTN(psimdata->pold), pbar = p - STARTP(psimdata->pold);

  double c   = PROCESS_GETVALUE_INSIDE(psimdata->c    , mbar, nbar, pbar);
  double rho = PROCESS_GETVALUE_INSIDE(psimdata->rho  , mbar, nbar, pbar);

  double dvx = PROCESS_GETVALUE_INSIDE(psimdata->vxold, mbar, nbar, pbar);
  double dvy = PROCESS_GETVALUE_INSIDE(psimdata->vyold, mbar, nbar, pbar);
  double dvz = PROCESS_GETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pbar);

  dvx -= m > 0 ? process_getvalue(psimdata->vxold, m - 1, n, p) : 0.0;
  dvy -= n > 0 ? process_getvalue(psimdata->vyold, m, n - 1, p) : 0.0;
  dvz -= p > 0 ? process_getvalue(psimdata->vzold, m, n, p - 1) : 0.0;

  double prev_p = PROCESS_GETVALUE_INSIDE(psimdata->pold, mbar, nbar, pbar);

  PROCESS_SETVALUE_INSIDE(psimdata->pnew, mbar, nbar, pbar, prev_p - rho * c * c * dtdx * (dvx + dvy + dvz));
}

void update_velocities(process_simulation_data_t *psimdata) {
  
  MPI_Request request_recv[3], request_send[3];
  register process_data_t *vxold = psimdata->vxold;

  int m, n, p;
  int mbar, nbar, pbar;
  register const int pnumnodesx = PNUMNODESX(vxold), startm = STARTM(vxold), endm = ENDM(vxold),
                     pnumnodesy = PNUMNODESY(vxold), startn = STARTN(vxold), endn = ENDN(vxold),
                     pnumnodesz = PNUMNODESZ(vxold), startp = STARTP(vxold), endp = ENDP(vxold);

  #pragma omp target teams distribute
  for (pbar = 0; pbar < pnumnodesz; pbar++) {
    #pragma omp parallel for
    for (nbar = 0; nbar < pnumnodesy; nbar++){
      psimdata->buffer_px[pbar * pnumnodesy + nbar] = PROCESS_GETVALUE_INSIDE(psimdata->pnew, 0, nbar, pbar);
    }
  }
  #pragma omp target teams distribute
  for (pbar = 0; pbar < pnumnodesz; pbar++) {    
    #pragma omp parallel for
    for (mbar = 0; mbar < pnumnodesx; mbar++){
      psimdata->buffer_py[pbar * pnumnodesx + mbar] = PROCESS_GETVALUE_INSIDE(psimdata->pnew, mbar, 0, pbar);
    }
  }
  #pragma omp target teams distribute
  for (nbar = 0; nbar < pnumnodesy; nbar++) {
    #pragma omp parallel for
    for (mbar = 0; mbar < pnumnodesx; mbar++){
      psimdata->buffer_pz[nbar * pnumnodesx + mbar] = PROCESS_GETVALUE_INSIDE(psimdata->pnew, mbar, nbar, 0);
    }
  }

  #pragma omp target update from(psimdata->buffer_px[0:pnumnodesy*pnumnodesz])
  #pragma omp target update from(psimdata->buffer_py[0:pnumnodesx*pnumnodesz])
  #pragma omp target update from(psimdata->buffer_pz[0:pnumnodesx*pnumnodesy])

  MPI_Isend(psimdata->buffer_px, pnumnodesy*pnumnodesz, MPI_DOUBLE, neighbors[LEFT ], SEND_X, cart_comm, &request_send[0]);
  MPI_Isend(psimdata->buffer_py, pnumnodesx*pnumnodesz, MPI_DOUBLE, neighbors[FRONT], SEND_Y, cart_comm, &request_send[1]);
  MPI_Isend(psimdata->buffer_pz, pnumnodesx*pnumnodesy, MPI_DOUBLE, neighbors[DOWN ], SEND_Z, cart_comm, &request_send[2]);

  MPI_Irecv(psimdata->pnew->ghostvals[RIGHT], pnumnodesy*pnumnodesz, MPI_DOUBLE, neighbors[RIGHT], SEND_X, cart_comm, &request_recv[0]);
  MPI_Irecv(psimdata->pnew->ghostvals[BACK ], pnumnodesx*pnumnodesz, MPI_DOUBLE, neighbors[BACK ], SEND_Y, cart_comm, &request_recv[1]);
  MPI_Irecv(psimdata->pnew->ghostvals[UP   ], pnumnodesx*pnumnodesy, MPI_DOUBLE, neighbors[UP   ], SEND_Z, cart_comm, &request_recv[2]);


  #pragma omp target teams distribute
  for (p = startp; p <= endp - 1; p++) {
    #pragma omp parallel for collapse(2)
    for (n = startn; n <= endn - 1; n++) {
      for (m = startm; m <= endm - 1; m++) {
        update_velocity_routine_inside(psimdata, m, n, p);
      }
    }
  }

  MPI_Waitall(3, request_recv, MPI_STATUS_IGNORE);

  #pragma omp target update to(psimdata->pnew->ghostvals[RIGHT][0:pnumnodesy*pnumnodesz])
  #pragma omp target update to(psimdata->pnew->ghostvals[BACK ][0:pnumnodesx*pnumnodesz])
  #pragma omp target update to(psimdata->pnew->ghostvals[UP   ][0:pnumnodesx*pnumnodesy])

  #pragma omp target teams distribute
  for (p = startp; p <= endp; p++) {
    #pragma omp parallel for
    for (n = startn; n <= endn; n++) {
      update_velocity_routine(psimdata, endm, n, p);
    }
  }
  #pragma omp target teams distribute
  for (p = startp; p <= endp; p++) {
    #pragma omp parallel for
    for (m = startm; m < endm; m++) {
      update_velocity_routine(psimdata, m, endn, p);
    }
  }
  #pragma omp target teams distribute
  for (n = startn; n < endn; n++){
    #pragma omp parallel for
    for (m = startm; m < endm; m++){
      update_velocity_routine(psimdata, m, n, endp);
    }
  }

  MPI_Waitall(3, request_send, MPI_STATUS_IGNORE);
}

void update_velocity_routine_inside(process_simulation_data_t *psimdata, int m, int n, int p) {

  process_data_t *pnew = psimdata->pnew;
  const register int mbar = m - STARTM(pnew), nbar = n - STARTN(pnew), pbar = p - STARTP(pnew);

  const double dtdxrho = psimdata->params.dt / psimdata->params.dx / PROCESS_GETVALUE_INSIDE(psimdata->rhohalf, mbar, nbar, pbar);

  double p_mnp = PROCESS_GETVALUE_INSIDE(pnew, mbar, nbar, pbar);

  double dpx = PROCESS_GETVALUE_INSIDE(pnew, mbar+1, nbar, pbar) - p_mnp;
  double dpy = PROCESS_GETVALUE_INSIDE(pnew, mbar, nbar+1, pbar) - p_mnp;
  double dpz = PROCESS_GETVALUE_INSIDE(pnew, mbar, nbar, pbar+1) - p_mnp;

  double prev_vx = PROCESS_GETVALUE_INSIDE(psimdata->vxold, mbar, nbar, pbar);
  double prev_vy = PROCESS_GETVALUE_INSIDE(psimdata->vyold, mbar, nbar, pbar);
  double prev_vz = PROCESS_GETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pbar);

  PROCESS_SETVALUE_INSIDE(psimdata->vxnew, mbar, nbar, pbar, prev_vx - dtdxrho * dpx);
  PROCESS_SETVALUE_INSIDE(psimdata->vynew, mbar, nbar, pbar, prev_vy - dtdxrho * dpy);
  PROCESS_SETVALUE_INSIDE(psimdata->vznew, mbar, nbar, pbar, prev_vz - dtdxrho * dpz);
}

void update_velocity_routine(process_simulation_data_t *psimdata, int m, int n, int p) {
  
  process_data_t *pnew = psimdata->pnew;
  const register int mbar = m - STARTM(pnew), nbar = n - STARTN(pnew), pbar = p - STARTP(pnew);

  const double dtdxrho = psimdata->params.dt / psimdata->params.dx / PROCESS_GETVALUE_INSIDE(psimdata->rhohalf, mbar, nbar, pbar);

  int mp1 = MIN(psimdata->global_grid.numnodesx - 1, m + 1);
  int np1 = MIN(psimdata->global_grid.numnodesy - 1, n + 1);
  int pp1 = MIN(psimdata->global_grid.numnodesz - 1, p + 1);

  double p_mnp = PROCESS_GETVALUE_INSIDE(pnew, mbar, nbar, pbar);

  double dpx = process_getvalue(pnew, mp1, n, p) - p_mnp;
  double dpy = process_getvalue(pnew, m, np1, p) - p_mnp;
  double dpz = process_getvalue(pnew, m, n, pp1) - p_mnp;

  double prev_vx = PROCESS_GETVALUE_INSIDE(psimdata->vxold, mbar, nbar, pbar);
  double prev_vy = PROCESS_GETVALUE_INSIDE(psimdata->vyold, mbar, nbar, pbar);
  double prev_vz = PROCESS_GETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pbar);

  PROCESS_SETVALUE_INSIDE(psimdata->vxnew, mbar, nbar, pbar, prev_vx - dtdxrho * dpx);
  PROCESS_SETVALUE_INSIDE(psimdata->vynew, mbar, nbar, pbar, prev_vy - dtdxrho * dpy);
  PROCESS_SETVALUE_INSIDE(psimdata->vznew, mbar, nbar, pbar, prev_vz - dtdxrho * dpz);
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
    MPI_Abort(cart_comm, MPI_ERR_IO);
  }

  if (cfp == NULL || c_numstep <= 0) {
    printf("Failed to open the speed map file. Aborting...\n\n");
    MPI_Abort(cart_comm, MPI_ERR_IO);
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
    MPI_Abort(cart_comm, MPI_ERR_IO);
  }

  fclose(rhofp);
  fclose(cfp);

  double xmin = rhoin_grid.xmin, xmax = rhoin_grid.xmax, 
    ymin = rhoin_grid.ymin, ymax = rhoin_grid.ymax, 
    zmin = rhoin_grid.zmin, zmax = rhoin_grid.zmax;

  double global_dx = xmax - xmin,
    global_dy = ymax - ymin,
    global_dz = zmax - zmin;

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

  psim_grid.lm.start = numnodesx*coords[0]/dims[0]; psim_grid.lm.end = numnodesx*(coords[0]+1)/dims[0] - 1;
  psim_grid.lm.n = psim_grid.lm.end - psim_grid.lm.start + 1;
  psim_grid.xmin = xmin + global_dx * psim_grid.lm.start / numnodesx;
  psim_grid.xmax = xmin + global_dx * (psim_grid.lm.end+1) / numnodesx;

  psim_grid.ln.start = numnodesy*coords[1]/dims[1]; psim_grid.ln.end = numnodesy*(coords[1]+1)/dims[1] - 1;
  psim_grid.ln.n = psim_grid.ln.end - psim_grid.ln.start + 1;
  psim_grid.ymin = ymin + global_dy * psim_grid.ln.start / numnodesy;
  psim_grid.ymax = ymin + global_dy * (psim_grid.ln.end+1) / numnodesy;

  psim_grid.lp.start = numnodesz*coords[2]/dims[2]; psim_grid.lp.end = numnodesz*(coords[2]+1)/dims[2] - 1;
  psim_grid.lp.n = psim_grid.lp.end - psim_grid.lp.start + 1;
  psim_grid.zmin = zmin + global_dz * psim_grid.lp.start / numnodesz;
  psim_grid.zmax = zmin + global_dz * (psim_grid.lp.end+1) / numnodesz;

  if (interpolate_inputmaps(psimdata, &psim_grid, c_map, rho_map) != 0) {
    printf(
        "Error while converting input map to simulation grid. Aborting...\n\n");
    MPI_Abort(cart_comm, MPI_ERR_NO_MEM);
  }

  DEBUG_PRINTF("Rank %4d has subdomain (%3d, %3d) (%3d, %3d) (%3d, %3d) of global grid (%3d, %3d, %3d)\n\
    \tsubdomain coords (%5.2lf, %5.2lf) (%5.2lf, %5.2lf) (%5.2lf, %5.2lf)\n\
    \tcvalue %10.5lf at (%3d, %3d, %3d)",
    cart_rank, psim_grid.lm.start, psim_grid.lm.end, psim_grid.ln.start, psim_grid.ln.end, psim_grid.lp.start, psim_grid.lp.end,
    psim_grid.gnumx, psim_grid.gnumy, psim_grid.gnumz, 
    psim_grid.xmin, psim_grid.xmax, psim_grid.ymin, psim_grid.ymax, psim_grid.zmin, psim_grid.zmax,
    psimdata->c->vals[0], psim_grid.lm.start, psim_grid.ln.start, psim_grid.lp.start);

  if (psimdata->params.outrate > 0 && psimdata->params.outputs != NULL) {
    for (int i = 0; i < psimdata->params.numoutputs; i++) {
      char *outfilei = psimdata->params.outputs[i].filename;

      for (int j = 0; j < i; j++) {
        char *outfilej = psimdata->params.outputs[j].filename;

        if (strcmp(outfilei, outfilej) == 0) {
          printf("Duplicate output file: '%s'. Aborting...\n\n", outfilei);
          MPI_Abort(cart_comm, MPI_ERR_IO);
        }
      }
    }

    for (int i = 0; i < psimdata->params.numoutputs; i++) {
      output_t *output = &psimdata->params.outputs[i];
      int m, n, p;
      closest_index(&psimdata->global_grid, output->posx, output->posy, output->posz, &m, &n, &p);

      if (open_outputfile(output, &psim_grid, &psimdata->global_grid) != 0) {
        printf("Failed to open output file: '%s'. Aborting...\n\n",
              output->filename);
        MPI_Abort(cart_comm, MPI_ERR_IO);
        break;
      }
    }
  }

  if ((psimdata->pold  = allocate_pdata(&psim_grid, 0b101010)) == NULL ||
      (psimdata->pnew  = allocate_pdata(&psim_grid, 0b101010)) == NULL ||
      (psimdata->vxold = allocate_pdata(&psim_grid, 0b000001)) == NULL ||
      (psimdata->vxnew = allocate_pdata(&psim_grid, 0b000001)) == NULL ||
      (psimdata->vyold = allocate_pdata(&psim_grid, 0b000100)) == NULL ||
      (psimdata->vynew = allocate_pdata(&psim_grid, 0b000100)) == NULL ||
      (psimdata->vzold = allocate_pdata(&psim_grid, 0b010000)) == NULL ||
      (psimdata->vznew = allocate_pdata(&psim_grid, 0b010000)) == NULL ||
      (psimdata->buffer_vx = malloc(psim_grid.ln.n * psim_grid.lp.n * sizeof(double))) == NULL ||
      (psimdata->buffer_vy = malloc(psim_grid.lm.n * psim_grid.lp.n * sizeof(double))) == NULL ||
      (psimdata->buffer_vz = malloc(psim_grid.lm.n * psim_grid.ln.n * sizeof(double))) == NULL ||
      (psimdata->buffer_px = malloc(psim_grid.ln.n * psim_grid.lp.n * sizeof(double))) == NULL ||
      (psimdata->buffer_py = malloc(psim_grid.lm.n * psim_grid.lp.n * sizeof(double))) == NULL ||
      (psimdata->buffer_pz = malloc(psim_grid.lm.n * psim_grid.ln.n * sizeof(double))) == NULL) {
    printf("Failed to allocate memory. Aborting...\n\n");
    MPI_Abort(cart_comm, MPI_ERR_NO_MEM);
  }

  fill_data(psimdata->pold, 0.0);
  fill_data(psimdata->pnew, 0.0);

  fill_data(psimdata->vynew, 0.0);
  fill_data(psimdata->vxold, 0.0);
  fill_data(psimdata->vynew, 0.0);
  fill_data(psimdata->vyold, 0.0);
  fill_data(psimdata->vznew, 0.0);
  fill_data(psimdata->vzold, 0.0);

  for(int i=0; i < PNUMNODESY(psimdata->pnew)*PNUMNODESZ(psimdata->pnew); i++){
    psimdata->pnew->ghostvals[RIGHT][i] = 0;
  }

  if(cart_rank == 0){
    printf("\n");
    printf(" Grid spacing: %g\n", psimdata->params.dx);
    printf("  Grid size X: %d\n", numnodesx);
    printf("  Grid size Y: %d\n", numnodesy);
    printf("  Grid size Z: %d\n", numnodesz);
    printf("    Time step: %g\n", psimdata->params.dt);
    printf(" Maximum time: %g\n\n", psimdata->params.maxt);

    if (psimdata->params.outrate > 0 && psimdata->params.outputs) {
      int outsampling =
          (int)(1.0 / (psimdata->params.outrate * psimdata->params.dt));

      printf("     Output rate: every %d step(s)\n", psimdata->params.outrate);
      printf(" Output sampling: %d Hz\n\n", outsampling);
      printf(" Output files:\n\n");

      for (int i = 0; i < psimdata->params.numoutputs; i++) {
        print_output(&psimdata->params.outputs[i]);
      }

      printf("\n");

    } else if (psimdata->params.outrate < 0) {
      printf("  Output is disabled (output rate set to 0)\n\n");

    } else {
      printf("  Output is disabled (not output specified)\n\n");
    }

    print_source(&psimdata->params.source);
    
    fflush(stdout);
  }

  free(rho_map->vals);
  free(rho_map);
  free(c_map->vals);
  free(c_map);
}

void finalize_simulation(process_simulation_data_t *psimdata) {
  if (psimdata->params.outputs != NULL) {
    for (int i = 0; i < psimdata->params.numoutputs; i++) {
      free(psimdata->params.outputs[i].filename);
      
      if (psimdata->params.outrate > 0 && psimdata->params.outputs[i].fp != NULL) {
        fclose(psimdata->params.outputs[i].fp);
      }
    }

    free(psimdata->params.outputs);
  }

  free(psimdata->params.source.data);
  free(psimdata->params.cin_filename);
  free(psimdata->params.rhoin_filename);

  free_pdata(psimdata->rho);
  free_pdata(psimdata->rhohalf);
  free_pdata(psimdata->c);

  free_pdata(psimdata->pold);
  free_pdata(psimdata->pnew);
  free_pdata(psimdata->vxold);
  free_pdata(psimdata->vyold);
  free_pdata(psimdata->vzold);
  free_pdata(psimdata->vxnew);
  free_pdata(psimdata->vynew);
  free_pdata(psimdata->vznew);

  free(psimdata->buffer_vx);
  free(psimdata->buffer_vy);
  free(psimdata->buffer_vz);
  free(psimdata->buffer_px);
  free(psimdata->buffer_py);
  free(psimdata->buffer_pz);
}

void swap_timesteps(process_simulation_data_t *psimdata) {

  const int numnodesx = PNUMNODESX(psimdata->pold);
  const int numnodesy = PNUMNODESY(psimdata->pold);
  const int numnodesz = PNUMNODESZ(psimdata->pold);

  #pragma omp target teams distribute
  for(int pbar = 0; pbar < numnodesz; pbar++){
    #pragma omp parallel for collapse(2)
    for(int nbar = 0; nbar < numnodesy; nbar++){
      for(int mbar = 0; mbar < numnodesx; mbar++){
        PROCESS_SETVALUE_INSIDE(psimdata->pold, mbar, nbar, pbar, 
          PROCESS_GETVALUE_INSIDE(psimdata->pnew, mbar, nbar, pbar));
        PROCESS_SETVALUE_INSIDE(psimdata->vxold, mbar, nbar, pbar, 
          PROCESS_GETVALUE_INSIDE(psimdata->vxnew, mbar, nbar, pbar));
        PROCESS_SETVALUE_INSIDE(psimdata->vyold, mbar, nbar, pbar, 
          PROCESS_GETVALUE_INSIDE(psimdata->vynew, mbar, nbar, pbar));
        PROCESS_SETVALUE_INSIDE(psimdata->vzold, mbar, nbar, pbar, 
          PROCESS_GETVALUE_INSIDE(psimdata->vznew, mbar, nbar, pbar));
      }
    }
  }
}
