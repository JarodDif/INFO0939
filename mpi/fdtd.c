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

  //Set correct value ranges
  // We have rankX, rankY, rankZ, numX, numY, numZ, *min/max
  psim_grid.gnumx = numnodesx; psim_grid.gnumy = numnodesy; psim_grid.gnumz = numnodesz;
  psim_grid.startm = numnodesx*coords[0]/dims[0]; psim_grid.endm = numnodesx*(coords[0]+1)/dims[0] - 1;
  psim_grid.startn = numnodesy*coords[1]/dims[1]; psim_grid.endn = numnodesy*(coords[1]+1)/dims[1] - 1;
  psim_grid.startp = numnodesz*coords[2]/dims[2]; psim_grid.endp = numnodesz*(coords[2]+1)/dims[2] - 1;

  printf("Rank %4d has has subdomain (%3d, %3d) (%3d, %3d) (%3d, %3d) of global grid (%3d, %3d, %3d)\n",
    cart_rank, psim_grid.startm, psim_grid.endm, psim_grid.startn, psim_grid.endn, psim_grid.startp, psim_grid.endp,
    psim_grid.gnumx, psim_grid.gnumy, psim_grid.gnumz);

#if 0

  if (interpolate_inputmaps(simdata, &sim_grid, c_map, rho_map) != 0) {
    printf(
        "Error while converting input map to simulation grid. Aborting...\n\n");
    exit(1);
  }

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
#if 0
  free(simdata->rho->vals);
  free(simdata->rho);
  free(simdata->rhohalf->vals);
  free(simdata->rhohalf);
  free(simdata->c->vals);
  free(simdata->c);

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

double process_getvalue(process_data_t *pdata, int m, int n, int p) {
  //Verify we are inside the boundaries
  if(m < 0 || n < 0 || p < 0) 
    return 0.0;
  if(m >= pdata->grid.gnumx || n >= pdata->grid.gnumy || p >= pdata->grid.gnumz)
    return 0.0;
  
  int mbar = m - STARTM(pdata);
  int nbar = n - STARTN(pdata);
  int pbar = p - STARTP(pdata);

  //We ask for ghost cells from the LEFT neighbor

  if(mbar >= 0 && mbar <= PNUMNODESX(pdata) - 1 &&
    nbar >= 0 && nbar <= PNUMNODESY(pdata) - 1 && 
    pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
    return pdata->vals[pbar*PNUMNODESX(pdata)*PNUMNODESY(pdata) + nbar*PNUMNODESX(pdata) + mbar];
  }

  if(mbar == - 1){
    if(nbar >= 0 && nbar <= PNUMNODESY(pdata) - 1 && pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return pdata->ghostvals[LEFT][pbar*PNUMNODESY(pdata) + nbar];
    }
    else{ return 0.0; }
  }else if(mbar == PNUMNODESX(pdata)){
    if(nbar >= 0 && nbar <= PNUMNODESY(pdata) - 1 && pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return pdata->ghostvals[RIGHT][pbar*PNUMNODESY(pdata) + nbar];
    }
    else{ return 0.0; }
  }else if(nbar == -1){
    if(pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return pdata->ghostvals[FRONT][pbar*PNUMNODESX(pdata) + mbar];
    }
    else{ return 0.0; }
  }else if(nbar == PNUMNODESY(pdata)){
    if(pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return pdata->ghostvals[BACK][pbar*PNUMNODESX(pdata) + mbar];
    }
    else{ return 0.0; }
  }else if(pbar == -1){
    return pdata->ghostvals[DOWN][nbar*PNUMNODESX(pdata) + mbar];
  }else if(pbar == PNUMNODESZ(pdata)){
    return pdata->ghostvals[UP][nbar*PNUMNODESX(pdata) + mbar];
  }
}

data_t *allocate_data(grid_t *grid) {
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

/******************************************************************************
 * Parameter file functions                                                   *
 ******************************************************************************/

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