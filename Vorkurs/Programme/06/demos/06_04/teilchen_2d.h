#ifndef TEILCHEN_2D_H
#define TEILCHEN_2D_H

typedef struct teilchen_2d_t {
  double x;
  double y;
  double v_x;
  double v_y;
  double m;
  double ladung;
} teilchen_2d_t;

void teilchen_2d_init(teilchen_2d_t * teilchen,
                      const unsigned int N);

void teilchen_2d_print(const teilchen_2d_t * teilchen);

#endif // ifndef(TEILCHEN_2D_H)
