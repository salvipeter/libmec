#pragma once

/* points is 2n long */
struct mec_t *mec_init(int n, const double *points);

/* point is 2 long, coordinates is n long */
void mec_eval(struct mec_t *mec, const double *point, double *coordinates);

void mec_free(struct mec_t *mec);
