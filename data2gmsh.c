#include <stdlib.h>
#include <string.h>

#include "datareader.h"

int write_data_gmsh(data_t *data, const char *filename, int step) {

  double delta_x =
      (NUMNODESX(data) > 1) ? ((XMAX(data) - XMIN(data)) / NUMNODESX(data)) : 0;
  double delta_y =
      (NUMNODESY(data) > 1) ? ((YMAX(data) - YMIN(data)) / NUMNODESY(data)) : 0;
  double delta_z =
      (NUMNODESZ(data) > 1) ? ((ZMAX(data) - ZMIN(data)) / NUMNODESZ(data)) : 0;

  FILE *fp;
  if ((fp = fopen(filename, step == 0 ? "w" : "a")) == NULL) {
    DEBUG_PRINTF("Could not open Gmsh data file '%s'", filename);
    return 1;
  }

  if (step == 0) {
    fprintf(fp, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    fprintf(fp, "$Nodes\n%d\n", NUMNODESTOT(data->grid));

    for (int p = 0; p < NUMNODESZ(data); p++) {
      for (int n = 0; n < NUMNODESY(data); n++) {
        for (int m = 0; m < NUMNODESX(data); m++) {
          int t =
              NUMNODESY(data) * NUMNODESX(data) * p + NUMNODESX(data) * n + m;
          double x = XMIN(data) + m * delta_x;
          double y = YMIN(data) + n * delta_y;
          double z = ZMIN(data) + p * delta_z;

          fprintf(fp, "%d %g %g %g\n", t, x, y, z);
        }
      }
    }

    fprintf(fp, "$EndNodes\n$Elements\n");

    if (NUMNODESX(data) == 1) {
      fprintf(fp, "%d\n", (NUMNODESY(data) - 1) * (NUMNODESZ(data) - 1));
      int e = 1;
      for (int p = 0; p < NUMNODESZ(data) - 1; p++) {
        for (int n = 0; n < NUMNODESY(data) - 1; n++) {
          int t1 = NUMNODESY(data) * NUMNODESX(data) * p + NUMNODESX(data) * n;
          int t2 =
              NUMNODESY(data) * NUMNODESX(data) * (p + 1) + NUMNODESX(data) * n;
          int t3 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) +
                   NUMNODESX(data) * (n + 1);
          int t4 =
              NUMNODESY(data) * NUMNODESX(data) * p + NUMNODESX(data) * (n + 1);

          fprintf(fp, "%d 3 2 1 1 %d %d %d %d\n", e++, t1, t2, t3, t4);
        }
      }
    } else if (NUMNODESY(data) == 1) {
      fprintf(fp, "%d\n", (NUMNODESX(data) - 1) * (NUMNODESZ(data) - 1));
      int e = 1;
      for (int p = 0; p < NUMNODESZ(data) - 1; p++) {
        for (int m = 0; m < NUMNODESX(data) - 1; m++) {
          int t1 = NUMNODESY(data) * NUMNODESX(data) * p + m;
          int t2 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) + m;
          int t3 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) + (m + 1);
          int t4 = NUMNODESY(data) * NUMNODESX(data) * p + (m + 1);

          fprintf(fp, "%d 3 2 1 1 %d %d %d %d\n", e++, t1, t2, t3, t4);
        }
      }
    } else if (NUMNODESZ(data) == 1) {
      fprintf(fp, "%d\n", (NUMNODESX(data) - 1) * (NUMNODESY(data) - 1));
      int e = 1;
      for (int n = 0; n < NUMNODESY(data) - 1; n++) {
        for (int m = 0; m < NUMNODESX(data) - 1; m++) {
          int t1 = NUMNODESX(data) * n + m;
          int t2 = NUMNODESX(data) * (n + 1) + m;
          int t3 = NUMNODESX(data) * (n + 1) + (m + 1);
          int t4 = NUMNODESX(data) * n + (m + 1);

          fprintf(fp, "%d 3 2 1 1 %d %d %d %d\n", e++, t1, t2, t3, t4);
        }
      }
    } else {
      fprintf(fp, "%d\n",
              (NUMNODESX(data) - 1) * (NUMNODESZ(data) - 1) *
                  (NUMNODESY(data) - 1));
      int e = 1;
      for (int p = 0; p < NUMNODESZ(data) - 1; p++) {
        for (int n = 0; n < NUMNODESY(data) - 1; n++) {
          for (int m = 0; m < NUMNODESX(data) - 1; m++) {
            int t1 =
                NUMNODESY(data) * NUMNODESX(data) * p + NUMNODESX(data) * n + m;
            int t2 = NUMNODESY(data) * NUMNODESX(data) * p +
                     NUMNODESX(data) * (n + 1) + m;
            int t3 = NUMNODESY(data) * NUMNODESX(data) * p +
                     NUMNODESX(data) * (n + 1) + (m + 1);
            int t4 = NUMNODESY(data) * NUMNODESX(data) * p +
                     NUMNODESX(data) * n + (m + 1);
            int t5 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) +
                     NUMNODESX(data) * n + m;
            int t6 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) +
                     NUMNODESX(data) * (n + 1) + m;
            int t7 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) +
                     NUMNODESX(data) * (n + 1) + (m + 1);
            int t8 = NUMNODESY(data) * NUMNODESX(data) * (p + 1) +
                     NUMNODESX(data) * n + (m + 1);

            fprintf(fp, "%d 5 2 1 1 %d %d %d %d %d %d %d %d\n", e++, t1, t2, t3,
                    t4, t5, t6, t7, t8);
          }
        }
      }
    }
    fprintf(fp, "$EndElements\n");
  }

  fprintf(fp, "$NodeData\n1\n\"%s\"\n1\n%d\n3\n%d\n1\n%d\n", filename, step,
          step, NUMNODESTOT(data->grid));

  for (int p = 0; p < NUMNODESZ(data); p++) {
    for (int n = 0; n < NUMNODESY(data); n++) {
      for (int m = 0; m < NUMNODESX(data); m++) {
        int t = NUMNODESY(data) * NUMNODESX(data) * p + NUMNODESX(data) * n + m;
        fprintf(fp, "%d %g\n", t, GETVALUE(data, m, n, p));
      }
    }
  }

  fprintf(fp, "$EndNodeData\n");
  fclose(fp);

  return 0;
}

data_t *extract_data(int operation, double xyz[3], data_t *in) {
  int m = 0, n = 0, p = 0;
  closest_index(in, xyz[0], xyz[1], xyz[2], &m, &n, &p);

  grid_descriptor_t grid;

  if (operation > 0 && operation < 4) {
    grid.numnodesx = (operation == 1) ? 1 : NUMNODESX(in);
    grid.numnodesy = (operation == 2) ? 1 : NUMNODESY(in);
    grid.numnodesz = (operation == 3) ? 1 : NUMNODESZ(in);

    grid.xmin = (operation == 1) ? xyz[0] : XMIN(in);
    grid.xmax = (operation == 1) ? xyz[0] : XMAX(in);
    grid.ymin = (operation == 2) ? xyz[1] : YMIN(in);
    grid.ymax = (operation == 2) ? xyz[1] : YMAX(in);
    grid.zmin = (operation == 3) ? xyz[2] : ZMIN(in);
    grid.zmax = (operation == 3) ? xyz[2] : ZMAX(in);

  } else {
    DEBUG_PRINTF("Unknown or invalid data extraction operation '%d'",
                 operation);
    return NULL;
  }

  data_t *data;
  if ((data = allocate_data(grid)) == NULL) {
    return NULL;
  }

  if (operation == 1) {
    for (p = 0; p < NUMNODESZ(in); p++) {
      for (n = 0; n < NUMNODESY(in); n++) {
        SETVALUE(data, 0, n, p, GETVALUE(in, m, n, p));
      }
    }

  } else if (operation == 2) {
    for (p = 0; p < NUMNODESZ(in); p++) {
      for (m = 0; m < NUMNODESX(in); m++) {
        SETVALUE(data, m, 0, p, GETVALUE(in, m, n, p));
      }
    }

  } else if (operation == 3) {
    for (n = 0; n < NUMNODESY(in); n++) {
      for (m = 0; m < NUMNODESX(in); m++) {
        SETVALUE(data, m, n, 0, GETVALUE(in, m, n, p));
      }
    }
  }

  return data;
}

void print_help(char *execname) {
  printf("Usage: %s [option] file(s)\n", execname);
  printf("Option can be\n");
  printf("  -cutx valx: to extract data on plane x = valx\n");
  printf("  -cuty valy: to extract data on plane y = valy\n");
  printf("  -cutz valz: to extract data on plane z = valz\n");
  printf("If no option is given, extract the full data\n\n");
  printf("The output .msh file can be opened with Gmsh (https://gmsh.info)\n");
}

int main(int argc, char **argv) {
  if (argc < 2) {
    print_help(argv[0]);
    return 1;
  }

  int operation = 0;
  double xyz[3] = {0., 0., 0.};

  char *filename = NULL;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-cutx") == 0 && i < argc - 1) {
      operation = 1;
      xyz[0] = atof(argv[++i]);

    } else if (strcmp(argv[i], "-cuty") == 0 && i < argc - 1) {
      operation = 2;
      xyz[1] = atof(argv[++i]);

    } else if (strcmp(argv[i], "-cutz") == 0 && i < argc - 1) {
      operation = 3;
      xyz[2] = atof(argv[++i]);

    } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      print_help(argv[0]);

      return 0;

    } else if (i == argc - 1) {
      filename = argv[i];

    } else {
      printf("Unknown option: '%s'. Aborting...\n\n", argv[i]);
      exit(1);
    }
  }

  int numsteps;
  grid_descriptor_t ingrid;

  FILE *fp = open_datafile(&ingrid, &numsteps, filename);

  if (fp == NULL) {
    printf("Failed to open data file (doesn't exist or is corrupted). Aborting...\n\n");
    exit(1);
  }

  if (operation == 0) {
    for (int step = 0; step < numsteps; step++) {
      data_t *data;

      if ((data = read_data(fp, ingrid, NULL, NULL)) == NULL) {
        printf("Failed to read data file. Aborting...\n\n");
        exit(1);
      }

      if (write_data_gmsh(data, "full.msh", step) != 0) {
        printf("Failed to write Gmsh data. Aborting...\n\n");
        exit(1);
      }

      free(data->vals);
      free(data);
    }

  } else {
    const char *outfile = (operation == 1)   ? "cutx.msh"
                          : (operation == 2) ? "cuty.msh"
                                             : "cutz.msh";

    for (int step = 0; step < numsteps; step++) {
      data_t *data, *cutdata;

      if ((data = read_data(fp, ingrid, NULL, NULL)) == NULL) {
        printf("Failed to read data file. Aborting...\n\n");
        exit(1);
      }

      if ((cutdata = extract_data(operation, xyz, data)) == NULL) {
        printf("Failed to extract data. Aborting...\n\n");
        exit(1);
      }

      free(data->vals);
      free(data);

      if (write_data_gmsh(cutdata, outfile, step) != 0) {
        printf("Failed to write Gmsh data. Aborting...\n\n");
        exit(1);
      }

      free(cutdata->vals);
      free(cutdata);
    }
  }

  printf("Done!\n");

  return 0;
}
