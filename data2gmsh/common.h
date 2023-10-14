#ifndef COMMON_H
#define COMMON_H

#if DEBUG

#define DEBUG_PRINTF(fmt, ...) printf("[DEBUG][%s:%d] " fmt "\n", __FILE__, __LINE__, __VA_ARGS__)
#define DEBUG_PRINT(msg)       printf("[DEBUG][%s:%d] %s\n", __FILE__, __LINE__, msg)

#else

#define DEBUG_PRINTF(fmt, ...)
#define DEBUG_PRINT(msg)

#define NUMNODESX(dat) (dat->grid.numnodesx)
#define NUMNODESY(dat) (dat->grid.numnodesy)
#define NUMNODESZ(dat) (dat->grid.numnodesz)

#define XMIN(dat) (dat->grid.xmin)
#define YMIN(dat) (dat->grid.ymin)
#define ZMIN(dat) (dat->grid.zmin)

#define XMAX(dat) (dat->grid.xmax)
#define YMAX(dat) (dat->grid.ymax)
#define ZMAX(dat) (dat->grid.zmax)

#define NUMNODESTOT(grid) (grid.numnodesx * grid.numnodesy * grid.numnodesz)

#define INDEX3D(grid, m, n, p) (grid.numnodesy * grid.numnodesx * (p) + \
                                                 grid.numnodesx * (n) + (m))

#define GETVALUE(dat, m, n, p)      (dat->vals[INDEX3D(dat->grid, m, n, p)])
#define SETVALUE(dat, m, n, p, val) (dat->vals[INDEX3D(dat->grid, m, n, p)] = val)

#endif

#endif