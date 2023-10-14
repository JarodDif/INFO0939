#include <mpi.h>
#include <stdio.h>

typedef enum neighbor {
  LEFT  = 0,
  RIGHT = 1,
  FRONT = 2,
  BACK  = 3,
  DOWN  = 4,
  UP    = 5,
} neighbor_t;

#define ndims 3

int main(int argc, char **argv) {
  int world_size;
  int rank, cart_rank;

  int dims[ndims] = {0, 0, 0};
  int periods[ndims] = {0, 0, 0};
  int reorder = 1;

  int coords[3];

  int neighbors[6];

  MPI_Comm cart_comm;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Dims_create(world_size, ndims, dims);

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm);
  MPI_Comm_rank(cart_comm, &cart_rank);

  MPI_Cart_coords(cart_comm, cart_rank, ndims, coords);

  MPI_Cart_shift(cart_comm, 0, 1, 
    &neighbors[LEFT], &neighbors[RIGHT]);

  MPI_Cart_shift(cart_comm, 1, 1, 
    &neighbors[FRONT], &neighbors[BACK]);

  MPI_Cart_shift(cart_comm, 2, 1, 
		&neighbors[DOWN], &neighbors[UP]);

  printf("Rank = %4d - Coords = (%3d, %3d, %3d) - Neighbors (up, down, left, right, front, back) = (%3d, %3d, %3d, %3d, %3d, %3d)\n",
    cart_rank, coords[0], coords[1], coords[2],
    neighbors[LEFT], neighbors[RIGHT], neighbors[FRONT], neighbors[BACK], neighbors[DOWN], neighbors[UP]);

  MPI_Finalize();

  return 0;
}

int main(int argc, const char *argv[]) {
  if (argc < 2) {
    printf("\nUsage: ./fdtd <param_file>\n\n");
    exit(1);
  }

  simulation_data_t simdata;
  init_simulation(&simdata, argv[1], coords);

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

  finalize_simulation(&simdata);

  return 0;
}

double process_getvalue(process_data_t *pdata, int m, int n, int p) {
  //Verify we are inside the boundaries
  if(m < 0 || n < 0 || p < 0) 
    return 0.0;
  if(m >= pdata.grid.gnumx || n >= pdata.grid.gnumy || p >= pdata.grid.gnumz)
    return 0.0;
  
  int mbar = m - STARTM(pdata);
  int nbar = n - STARTN(pdata);
  int pbar = p - STARTP(pdata);

  //We ask for ghost cells from the LEFT neighbor

  if(mbar >= 0 && mbar <= PNUMNODESX(pdata) - 1 &&
    nbar >= 0 && nbar <= PNUMNODESY(pdata) - 1 && 
    pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
    return pdata->vals[pbar*PNUMNODESX(pdata)*PNUMNODESY(pdata) + nbar*PNUMNODEX(pdata) + mbar];
  }

  if(mbar == - 1){
    if(nbar >= 0 && nbar <= PNUMNODESY(pdata) - 1 && pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return data->ghostvals[LEFT][pbar*PNUMNODESY(pdata) + nbar];
    }
    else{ return 0.0; }
  }else if(mbar == PNUMNODESX(pdata)){
    if(nbar >= 0 && nbar <= PNUMNODESY(pdata) - 1 && pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return data->ghostvals[RIGHT][pbar*PNUMNODESY(pdata) + nbar];
    }
    else{ return 0.0; }
  }else if(nbar == -1){
    if(pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return data->ghostvals[FRONT][pbar*PNUMNODEX(pdata) + mbar]
    }
    else{ return 0.0; }
  }else if(nbar == PNUMNODESY(pdata)){
    if(pbar >= 0 && pbar <= PNUMNODESZ(pdata) - 1){
      return data->ghostvals[BACK][pbar*PNUMNODESX(pdata) + mbar];
    }
    else{ return 0.0; }
  }else if(pbar == -1){
    return data->ghostvals[DOWN][nbar*PNUMNODEX(pdata) + mbar];
  }else if(pbar == PNUMNODESZ(pdata)){
    return data->ghostvals[UP][nbar*PNUMNODEX(pdata) + mbar];
  }
}

void init_simulation(simulation_data_t *simdata, const char *params_filename, const int* dims, const int* coords) {
  if (read_paramfile(&simdata->params, params_filename) != 0) {
    printf("Failed to read parameters. Aborting...\n\n");
    exit(1);
  }

  grid_t rhoin_grid;
  grid_t cin_grid;
  grid_t sim_grid;

  int rho_numstep;
  int c_numstep;

  FILE *rhofp =
      open_datafile(&rhoin_grid, &rho_numstep, simdata->params.rhoin_filename);
  FILE *cfp =
      open_datafile(&cin_grid, &c_numstep, simdata->params.cin_filename);

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

  /* TODO: We still need to set these values
  sim_grid.xmin = rhoin_grid.xmin;
  sim_grid.xmax = rhoin_grid.xmax;
  sim_grid.ymin = rhoin_grid.ymin;
  sim_grid.ymax = rhoin_grid.ymax;
  sim_grid.zmin = rhoin_grid.zmin;
  sim_grid.zmax = rhoin_grid.zmax;
  */

  int numnodesx = MAX(floor((xmax - xmin) / simdata->params.dx), 1),
    numnodesy = MAX(floor((ymax - ymin) / simdata->params.dx), 1),
    numnodesz = MAX(floor((zmax - zmin) / simdata->params.dx), 1);

  /* TODO: We need to set these values
  sim_grid.numnodesx =
      MAX(floor((sim_grid.xmax - sim_grid.xmin) / simdata->params.dx), 1);
  sim_grid.numnodesy =
      MAX(floor((sim_grid.ymax - sim_grid.ymin) / simdata->params.dx), 1);
  sim_grid.numnodesz =
      MAX(floor((sim_grid.zmax - sim_grid.zmin) / simdata->params.dx), 1);
  */

  //Set correct value ranges
  // We have rankX, rankY, rankZ, numX, numY, numZ, *min/max

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

  free(rho_map->vals);
  free(rho_map);
  free(c_map->vals);
  free(c_map);
}

void finalize_simulation(simulation_data_t *simdata) {
  if (simdata->params.outputs != NULL) {
    for (int i = 0; i < simdata->params.numoutputs; i++) {
      free(simdata->params.outputs[i].filename);

      if (simdata->params.outrate > 0) {
        fclose(simdata->params.outputs[i].fp);
      }
    }

    free(simdata->params.outputs);
  }

  free(simdata->params.source.data);
  free(simdata->params.cin_filename);
  free(simdata->params.rhoin_filename);

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
}

void swap_timesteps(simulation_data_t *simdata) {
  data_t *tmpp = simdata->pold;
  data_t *tmpvx = simdata->vxold;
  data_t *tmpvy = simdata->vyold;
  data_t *tmpvz = simdata->vzold;

  simdata->pold = simdata->pnew;
  simdata->pnew = tmpp;
  simdata->vxold = simdata->vxnew;
  simdata->vxnew = tmpvx;
  simdata->vyold = simdata->vynew;
  simdata->vynew = tmpvy;
  simdata->vzold = simdata->vznew;
  simdata->vznew = tmpvz;
}