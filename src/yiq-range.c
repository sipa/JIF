#include <stdio.h>
#include <stdlib.h>

#include "image/pixel.h"
#include "color.h"

char posy[256] = {};
char posi[256][511] = {};
char posq[256][511][511] = {};

int main(void) {
  for (int r = 0; r < 256; r++) {
    for (int b = 0; b < 256; b++) {
      for (int g = 0; g < 256; g++) {
        pixel_t p={.d={r,g,b}};
        pixel_rgb2yiq(&p);
        posy[UNSTRETCH(p.d[0])] = 1;
        posi[UNSTRETCH(p.d[0])][UNSTRETCH(p.d[1])] = 1;
        posq[UNSTRETCH(p.d[0])][UNSTRETCH(p.d[1])][UNSTRETCH(p.d[2])] = 1;
      }
    }
  }
  {
    int low = 255, high = 0;
    for (int y = 0; y < 256; y++) {
      if (posy[y]) {
        if (y < low) low = y;
        if (y > high) high = y;
      }
    }
    int min, max;
    get_range_y(&min,&max);
    min=UNSTRETCH(min);
    max=UNSTRETCH(max);
    if (min != low || max != high) {
      printf("Y=%i..%i (not %i..%i)\n",low,high,min,max);
    } else {
      //      printf("Y=%i..%i (ok)\n",low,high);
    }
  }
  for (int y = 0; y < 256; y++) {
    int low = 510, high = 0;
    for (int i = 0; i < 511; i++) {
      if (posi[y][i]) {
        if (i < low) low = i;
        if (i > high) high = i;
      }
    }
    int min, max;
    get_range_i(STRETCH(y),&min,&max);
    min=UNSTRETCH(min);
    max=UNSTRETCH(max);
    if (min != low || max != high) {
      printf("Y=%i: I=%i..%i (not %i..%i)\n",y,low,high,min,max);
    } else {
      //      printf("Y=%i: I=%i..%i (ok)\n",y,low,high);
    }
  }
  for (int y = 0; y < 256; y++) {
    for (int i = 0; i < 511; i++) {
      if (posi[y][i]) {
        int low = 510, high = 0;
        for (int q = 0; q < 511; q++) {
          if (posq[y][i][q]) {
            if (q < low) low = q;
            if (q > high) high = q;
          }
        }
        int min, max;
        get_range_q(STRETCH(y),STRETCH(i),&min,&max);
        min=UNSTRETCH(min);
        max=UNSTRETCH(max);
        if (min != low || max != high) {
          printf("Y=%i I=%i: Q=%i..%i (not %i..%i)\n",y,i,low,high,min,max);
        } else {
          //          printf("Y=%i I=%i: Q=%i..%i (ok)\n",y,i,low,high);
        }
      }
    }
  }
  return 0;
}
