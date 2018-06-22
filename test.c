#include <stdio.h>

#include "mec.h"

int main(int argc, char **argv) {
  double points[] = { -1, -1, 1, -1, 1, 1, 0.33, 1, 0.33, 0.33, -0.33, 0.33, -0.33, 1, -1, 1 };
  double point[] = { -1, -1 };
  double coordinates[8];
  struct mec_t *mec = mec_init(8, points);
  mec_eval(mec, point, coordinates);
  mec_free(mec);
  for (int i = 0; i < 8; ++i)
    printf("%lf\n", coordinates[i]);
  return 0;
}
