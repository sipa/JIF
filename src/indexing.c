// TODO: auto-"indexing" (moet nog verder geimplementeerd worden)

// idee is om veel nauwkeurigere get_range_{y,i,q}() te maken op basis van
// een soort indexing. Elke bucket komt overeen met een balkje in de Y(I(Q))-ruimte,
// en kan enkele specifieke kleuren bevatten. Als er teveel verschillende kleuren in
// 1 bucket vallen dan wordt de bucket "continue" en zijn alle kleuren daarin mogelijk.
// Als het aantal kleuren beperkt is, dan is de bucket "discreet" en zullen hopelijk veel
// ranges kunnen beperkt worden tot lege intervallen.

// de bedoeling is dat reduce_range() wordt gebruikt tijdens het outputten van pixel diffs
// om zo weinig mogelijk redundante bits te moeten outputten
// als reduce_range() een singleton interval geeft dan moet er verder niets worden geoutput
// als reduce_range() een leeg interval geeft dan is er geen kleur mogelijk, dus bits waarvan
// 1 vd 2 waarden een leeg interval veroorzaken zijn ook redundant


#include "config.h"
#include "util.h"
#include "image/pixel.h"
#include "indexing.h"
#include "bitvector.h"
#include <stdio.h>

// stretched a, unstretched b
int pixel_dist_su(const pixel_t *a, const pixel_t *b) {
  int ydiff = abs(UNSTRETCH(a->d[0])- b->d[0]);
  int idiff = abs(UNSTRETCH(a->d[1])- b->d[1]);
  int qdiff = abs(UNSTRETCH(a->d[2])- b->d[2]);
  return 2*ydiff + idiff + qdiff;
}

void round_color(const colors *a, pixel_t *pixel) {
  if (a->used==0) return;
  int y = UNSTRETCH(pixel->d[0]);
  int i = UNSTRETCH(pixel->d[1]);
  int q = UNSTRETCH(pixel->d[2]);
//  fprintf(stderr,"color rounding: (%i,%i,%i)\n",y,i,q);

  int kyb = clamp((y / SIZE_Y)-1,0,BUCKETS_Y-1);
  int kib = clamp((i / SIZE_I)-1,0,BUCKETS_I-1);
  int kqb = clamp((q / SIZE_Q)-1,0,BUCKETS_Q-1);
  int kye = clamp((y / SIZE_Y)+1,0,BUCKETS_Y-1);
  int kie = clamp((i / SIZE_I)+1,0,BUCKETS_I-1);
  int kqe = clamp((q / SIZE_Q)+1,0,BUCKETS_Q-1);

  int nearest = -1;
  pixel_t newpixelV ={};
  pixel_t * newpixel = &newpixelV;
  int no_unique_nearest = 0;

        // start with middle bucket, it may contain exact match
        int mky=y/SIZE_Y;
        int mki=i/SIZE_I;
        int mkq=q/SIZE_Q;
        int c = a->bucketYIQ[mky][mki][mkq].count;
        if (c > 0) {
          if (c < MAX_PER_BUCKET) {
            for (int j=0; j < c && nearest != 0; j++) {
                pixel_t * color = &a->bucketYIQ[mky][mki][mkq].color[j];
                int dist = pixel_dist_su(pixel, color);
                if (dist == nearest) no_unique_nearest++;
                if (nearest < 0 || dist < nearest) {
                        newpixel = color;
                        nearest = dist;
                        no_unique_nearest=0;
                }
            }
          } else {
            pixel_t best_color = {};
//            best_color.d[0] = clamp(y,ky*SIZE_Y,(ky+1)*SIZE_Y-1);
 //           best_color.d[1] = clamp(i,ki*SIZE_I,(ki+1)*SIZE_I-1);
   //         best_color.d[2] = clamp(q,kq*SIZE_Q,(kq+1)*SIZE_Q-1);
            best_color.d[0] = clamp(y,a->bucketYIQ[mky][mki][mkq].color[0].d[0],a->bucketYIQ[mky][mki][mkq].color[1].d[0]);
            best_color.d[1] = clamp(i,a->bucketYIQ[mky][mki][mkq].color[0].d[1],a->bucketYIQ[mky][mki][mkq].color[1].d[1]);
            best_color.d[2] = clamp(q,a->bucketYIQ[mky][mki][mkq].color[0].d[2],a->bucketYIQ[mky][mki][mkq].color[1].d[2]);
            int dist = pixel_dist_su(pixel, &best_color);
            if (dist == nearest) no_unique_nearest++;
            if (nearest < 0 || dist < nearest) {
                   newpixelV = best_color;
                   newpixel = &newpixelV;
                   nearest = dist;
                   no_unique_nearest=0;
            }
          }
         }


  for (int ky=kyb; ky <= kye && nearest != 0; ky++) {
    for (int ki=kib; ki <= kie && nearest != 0; ki++) {
      for (int kq=kqb; kq <= kqe && nearest != 0; kq++) {
        if (ky==mky && ki==mki && kq == mkq) continue;
        int c = a->bucketYIQ[ky][ki][kq].count;
        if (c > 0) {
          if (c < MAX_PER_BUCKET) {
            for (int j=0; j < c && nearest != 0; j++) {
                pixel_t * color = &a->bucketYIQ[ky][ki][kq].color[j];
                int dist = pixel_dist_su(pixel, color);
                if (dist == nearest) no_unique_nearest++;
                if (nearest < 0 || dist < nearest) {
                        newpixel = color;
                        nearest = dist;
                        no_unique_nearest=0;
                }
            }
          } else {
            pixel_t best_color = {};
//            best_color.d[0] = clamp(y,ky*SIZE_Y,(ky+1)*SIZE_Y-1);
 //           best_color.d[1] = clamp(i,ki*SIZE_I,(ki+1)*SIZE_I-1);
   //         best_color.d[2] = clamp(q,kq*SIZE_Q,(kq+1)*SIZE_Q-1);
            best_color.d[0] = clamp(y,a->bucketYIQ[ky][ki][kq].color[0].d[0],a->bucketYIQ[ky][ki][kq].color[1].d[0]);
            best_color.d[1] = clamp(i,a->bucketYIQ[ky][ki][kq].color[0].d[1],a->bucketYIQ[ky][ki][kq].color[1].d[1]);
            best_color.d[2] = clamp(q,a->bucketYIQ[ky][ki][kq].color[0].d[2],a->bucketYIQ[ky][ki][kq].color[1].d[2]);
            int dist = pixel_dist_su(pixel, &best_color);
            if (dist == nearest) no_unique_nearest++;
            if (nearest < 0 || dist < nearest) {
                   newpixelV = best_color;
                   newpixel = &newpixelV;
                   nearest = dist;
                   no_unique_nearest=0;
            }
          }
       }
      }
    }
  }
  if (nearest != -1 && no_unique_nearest==0) {
/*        if (nearest>0)
        fprintf(stderr,"color rounding: (%i,%i,%i) -> (%i,%i,%i)  [dist=%i]\n",
                UNSTRETCH(pixel->d[0]),
                UNSTRETCH(pixel->d[1]),
                UNSTRETCH(pixel->d[2]),
                newpixel->d[0],
                newpixel->d[1],
                newpixel->d[2],
                nearest);*/
        for(int c=0; c<3; c++) pixel->d[c] = STRETCH(newpixel->d[c]);
  }
}
int bucket_exists(int yk, int ik, int qk) {
  int ymin = yk*SIZE_Y;
  int ymax = (yk+1)*SIZE_Y-1;
  int imin = ik*SIZE_I;
  int imax = (ik+1)*SIZE_I-1;
  int cimin = 0;
  int cimax = 510;
  if (ymax<63) {
    cimin=252-4*ymax;
    cimax=258+4*ymax;
  } else if (ymin>=192) {
    cimin=3+4*(ymin-192);
    cimax=507-4*(ymin-192);
  }
  if (cimin > imax) return 0;
  if (cimax < imin) return 0;

  int qmin = qk*SIZE_Q;
  int qmax = (qk+1)*SIZE_Q-1;
  int cqmin = 0;
  int cqmax = 510;
  if (ymax<63) {
    cqmin=252-4*ymax;
    cqmax=258+4*ymax;
  } else if (ymin>=192) {
    cqmin=3+4*(ymin-192);
    cqmax=507-4*(ymin-192);
  }
  if (cqmin > qmax) return 0;
  if (cqmax < qmin) return 0;


/* berekening klopt niet
  int qmin = qk*SIZE_Q;
  int qmax = (qk+1)*SIZE_Q-1;
  int cqmin = 0;
  int cqmax = 510;
  int cimind = MIN(abs(cimin-255),abs(cimax-255));
//  int cimaxd = MAX(abs(cimin-255),abs(cimax-255));
  if (ymax<63) {
    cqmin=254-2*ymax+(cimind/2)*2;
    cqmax=256+2*ymax;
  } else if (ymin>=192) {
    cqmin=255-2*(255-ymin);
    cqmax=255+2*(255-ymin)-((1+cimind)/2)*2;
  } else {
//    cqmin=MIN(1+(ymin-128)*2,128-(ymax-63)*2+(cimind/2)*2);
//    cqmax=MAX(382+(ymax-63)*2,383+(191-ymin)*2-((1+cimind)/2)*2);
  }
  if (cqmin > qmax) return 0;
  if (cqmax < qmin) return 0;
*/
  return 1;
}


int similar(int x, int y, int d) {
  return x==y;
/*  if (x >= y+d) return 0;
  if (x <= y-d) return 0;
  return 1;*/
}
/*
int reduce_range_Y(colors *a, int *y_min, int *y_max) {
   int ky_min = *y_min / SIZE_Y;
   int ky_max = *y_max / SIZE_Y;
   int new_y_min = *y_max+1;
   int new_y_max = *y_min-1;
   for (int i=ky_min; i <= ky_max; i++) {
     assert(i<=BUCKETS_Y);
     int c = a->bucketY[i].count;
     if (c > 0) {
        if (c < MAX_PER_BUCKET) {
            for (int j=0; j < c; j++) {
                int y = a->bucketY[i].color[j].d[0];
                if (y >= *y_min && y <= *y_max) {
                    if (y < new_y_min) new_y_min = y;
                    if (y > new_y_max) new_y_max = y;
                }
            }
        } else {
//            new_y_min = clamp(i*SIZE_Y,*y_min,new_y_min);
//            new_y_max = clamp((i+1)*SIZE_Y-1,new_y_max,*y_max);
            assert(a->bucketY[i].color[0].d[0] >= i*SIZE_Y);
            assert(a->bucketY[i].color[1].d[0] <= (i+1)*SIZE_Y-1);
          new_y_min = clamp(a->bucketY[i].color[0].d[0],*y_min,new_y_min);
          new_y_max = clamp(a->bucketY[i].color[1].d[0],new_y_max,*y_max);
        }
     }
   }
   *y_min = new_y_min;
   *y_max = new_y_max;
   if (new_y_min == new_y_max) return 2;        // exact color known
   if (new_y_min < new_y_max) return 1;         // still several options
   return 0;                                    // no color found
}
int inline reduce_range_I(colors *a, int *i_min, int *i_max, int y) {
   int ki_min = *i_min / SIZE_I;
   int ki_max = *i_max / SIZE_I;
   int ky_min = (y-a->max_diff[0]) / SIZE_Y;
   int ky_max = (y+a->max_diff[0]) / SIZE_Y;
   int new_i_min = *i_max+1;
   int new_i_max = *i_min-1;
   for (int ky=ky_min; ky <= ky_max; ky++) {
    for (int i=ki_min; i <= ki_max; i++) {
     int c = a->bucketYI[ky][i].count;
     if (c > 0) {
        if (c < MAX_PER_BUCKET) {
            for (int j=0; j < c; j++) {
                int ty = a->bucketYI[ky][i].color[j].d[0];
                int ti = a->bucketYI[ky][i].color[j].d[1];
                if (similar(ty,y,a->max_diff[0]) && ti >= *i_min && ti <= *i_max) {
                    if (ti < new_i_min) new_i_min = ti;
                    if (ti > new_i_max) new_i_max = ti;
                }
            }
        } else {
//            new_i_min = clamp(i*SIZE_I,*i_min,new_i_min);
//            new_i_max = clamp((i+1)*SIZE_I-1,new_i_max,*i_max);
            assert(a->bucketYI[ky][i].color[0].d[1] >= i*SIZE_I);
            assert(a->bucketYI[ky][i].color[1].d[1] <= (i+1)*SIZE_I-1);
            int y_min = a->bucketYI[ky][i].color[0].d[0];
            int y_max = a->bucketYI[ky][i].color[1].d[0];
            if (y >= y_min && y <= y_max) {
              new_i_min = clamp(a->bucketYI[ky][i].color[0].d[1],*i_min,new_i_min);
              new_i_max = clamp(a->bucketYI[ky][i].color[1].d[1],new_i_max,*i_max);
            }
            // TODO: check alle bucketYIQ's voor betere grenzen
        }
     }
    }
   }
   *i_min = new_i_min;
   *i_max = new_i_max;
   if (new_i_min == new_i_max) return 2;        // exact color known
   if (new_i_min < new_i_max) return 1;         // still several options
   return 0;                                    // no color found
}
int inline reduce_range_Q(colors *a, int *q_min, int *q_max, int y, int mi) {
//   int kq_min = clamp(*q_min / SIZE_Q -1, 0, BUCKETS_Q-1);
//   int kq_max = clamp(*q_max / SIZE_Q + 1, 0, BUCKETS_Q-1);
   int kq_min = *q_min / SIZE_Q;
   int kq_max = *q_max / SIZE_Q;
   int ky_min = (y-a->max_diff[0]) / SIZE_Y;
   int ky_max = (y+a->max_diff[0]) / SIZE_Y;
   int ki_min = (mi-a->max_diff[1]) / SIZE_I;
   int ki_max = (mi+a->max_diff[1]) / SIZE_I;

   int new_q_min = *q_max+1;
   int new_q_max = *q_min-1;
   for (int ky=ky_min; ky <= ky_max; ky++) {
   for (int ki=ki_min; ki <= ki_max; ki++) {
    for (int i=kq_min; i <= kq_max; i++) {
     int c = a->bucketYIQ[ky][ki][i].count;
//     fprintf(stderr,"Looking in bucket[%i][%i][%i] (contains %i colors)\n",ky,ki,i,c);
     if (c > 0) {
        if (c < MAX_PER_BUCKET) {
            for (int j=0; j < c; j++) {
                int ty = a->bucketYIQ[ky][ki][i].color[j].d[0];
                int ti = a->bucketYIQ[ky][ki][i].color[j].d[1];
                int tq = a->bucketYIQ[ky][ki][i].color[j].d[2];
//                fprintf(stderr,"Considering color Y=%i,I=%i,Q=%i\n",ty,ti,tq);
                if (similar(ty,y,a->max_diff[0]) && similar(ti,mi,a->max_diff[1]) && tq >= *q_min && tq <= *q_max) {
                    if (tq < new_q_min) new_q_min = tq;
                    if (tq > new_q_max) new_q_max = tq;
                }
            }
        } else {
//            new_q_min = clamp(i*SIZE_Q,*q_min,new_q_min);
//            new_q_max = clamp((i+1)*SIZE_Q-1,new_q_max,*q_max);
            assert(a->bucketYIQ[ky][ki][i].color[0].d[2] >= i*SIZE_Q);
            assert(a->bucketYIQ[ky][ki][i].color[1].d[2] <= (i+1)*SIZE_Q-1);
            int y_min = a->bucketYIQ[ky][ki][i].color[0].d[0];
            int y_max = a->bucketYIQ[ky][ki][i].color[1].d[0];
            int i_min = a->bucketYIQ[ky][ki][i].color[0].d[1];
            int i_max = a->bucketYIQ[ky][ki][i].color[1].d[1];
            if (y >= y_min && y <= y_max && mi >= i_min && mi <= i_max) {
              new_q_min = clamp(a->bucketYIQ[ky][ki][i].color[0].d[2],*q_min,new_q_min);
              new_q_max = clamp(a->bucketYIQ[ky][ki][i].color[1].d[2],new_q_max,*q_max);
            }
//            if (*q_min != new_q_min || *q_max != new_q_max)
//            fprintf(stderr,"(y=%i,i=%i) old range: %i..%i, new range: %i..%i\n",y,ki,*q_min,*q_max,new_q_min,new_q_max);
        }
     }
    }
   }
   }
//   fprintf(stderr,"Y=%i, I=%i, Q in %i..%i -> %i..%i\n",y,mi,*q_min,*q_max,new_q_min,new_q_max);
   *q_min = new_q_min;
   *q_max = new_q_max;
   if (new_q_min == new_q_max) return 2;        // exact color known
   if (new_q_min < new_q_max) return 1;         // still several options
   return 0;                                    // no color found
}

int reduce_range_YIQ(colors *a, int *ranges, int c) {

   int ky_min = clamp(ranges[0] / SIZE_Y, 0, BUCKETS_Y-1);
   int ky_max = clamp(ranges[3] / SIZE_Y, 0, BUCKETS_Y-1);
   int ki_min = clamp(ranges[1] / SIZE_I, 0, BUCKETS_I-1);
   int ki_max = clamp(ranges[4] / SIZE_I, 0, BUCKETS_I-1);
   int kq_min = clamp(ranges[2] / SIZE_Q, 0, BUCKETS_Q-1);
   int kq_max = clamp(ranges[5] / SIZE_Q, 0, BUCKETS_Q-1);

   int new_min = ranges[3+c]+1;
   int new_max = ranges[c]-1;
   for (int ky=ky_min; ky <= ky_max; ky++) {
    if (a->bucketY[ky].count == 0) continue;
    for (int ki=ki_min; ki <= ki_max; ki++) {
     if (a->bucketYI[ky][ki].count == 0) continue;
     for (int kq=kq_min; kq <= kq_max; kq++) {
      int count = a->bucketYIQ[ky][ki][kq].count;
      if (count > 0) {
        if (count < MAX_PER_BUCKET) {
            for (int j=0; j < count; j++) {
                int inrange = 1;
                for (int i=0; i<3; i++) {
                  if (a->bucketYIQ[ky][ki][kq].color[j].d[i] < ranges[i]) inrange=0;
                  if (a->bucketYIQ[ky][ki][kq].color[j].d[i] > ranges[3+i]) inrange=0;
                }
                if (inrange) {
                    if (a->bucketYIQ[ky][ki][kq].color[j].d[c] < new_min) new_min = a->bucketYIQ[ky][ki][kq].color[j].d[c];
                    if (a->bucketYIQ[ky][ki][kq].color[j].d[c] > new_max) new_max = a->bucketYIQ[ky][ki][kq].color[j].d[c];
                }
            }
        } else {
            int inrange = 1;
            for (int i=0; i<3; i++) {
              if (a->bucketYIQ[ky][ki][kq].color[1].d[i] < ranges[i]) inrange=0;        // [bucket] [range]
              if (a->bucketYIQ[ky][ki][kq].color[0].d[i] > ranges[3+i]) inrange=0;      // [range] [bucket]
              // otherwise overlap
            }
            if (inrange) {
              new_min = clamp(a->bucketYIQ[ky][ki][kq].color[0].d[c],ranges[c],new_min);
              new_max = clamp(a->bucketYIQ[ky][ki][kq].color[1].d[c],new_max,ranges[3+c]);
            }
        }
     }
    }
   }
   }
   ranges[c] = new_min;
   ranges[3+c] = new_max;
   if (new_min == new_max) return 2;        // exact color known
   if (new_min < new_max) return 1;         // still several options
   return 0;                                // no color found
}

int inline reduce_range(colors *a, int *min, int *max, int c, pixel_t *known, int qz, int guess) {
  if (*min > *max) return 0;

  if (c<0 || a->used==0 || a->max_diff[c] > 1) return (*min < *max ? 1 : (*min == *max ? 2 : 0));
  int result = 1;
  int ranges[6] = {0,0,0,255,510,510};
  ranges[c] = clamp(UNSTRETCH(*min * qz + guess),0,(c==0?255:510));
  ranges[3+c] = clamp(UNSTRETCH(*max * qz + guess),0,(c==0?255:510));
  for(int i=0; i<c; i++) {
    ranges[i] = UNSTRETCH(known->d[i]);
    ranges[3+i] = UNSTRETCH(known->d[i]);
  }
//  fprintf(stderr,"reduced range (channel %i, guess %i, qz %i) from %i..%i (%i..%i) to ",c,guess,qz,*min,*max,ranges[c],ranges[3+c]);
  result = reduce_range_YIQ(a,&ranges[0],c);
  *min = (STRETCH(ranges[c])-guess)/qz;
  *max = (STRETCH(ranges[3+c])-guess)/qz;
//  fprintf(stderr,"%i..%i (%i..%i)  (result:%i)\n",*min,*max,ranges[c],ranges[3+c],result);
  return result;
}
*/




/*
int ranges_overlap_bucket(color_bucket *b, int ranges[6]) {
    int count = b->count;
    if (count > 0) {
         if (count < MAX_PER_BUCKET) {
             for (int j=0; j < count; j++) {
                 int inrange = 1;
                 for (int i=0; i<3; i++) {
                   if (b->color[j].d[i] < ranges[i]) inrange=0;
                   if (b->color[j].d[i] > ranges[3+i]) inrange=0;
                 }
                 if (inrange) return 1;
             }
         } else {
             int inrange = 1;
             for (int i=0; i<3; i++) {
               if (b->color[1].d[i] < ranges[i]) inrange=0;        // [bucket] [range]
               if (b->color[0].d[i] > ranges[3+i]) inrange=0;      // [range] [bucket]
               // otherwise overlap
             }
             if (inrange) return 1;
         }
    }
    return 0;
}



inline int color_exists(colors *a, int ranges[6], int channel) {
   int ky_min = ranges[0] / SIZE_Y;
   int ky_max = ranges[3] / SIZE_Y;
   int ki_min = ranges[1] / SIZE_I;
   int ki_max = ranges[4] / SIZE_I;
   int kq_min = ranges[2] / SIZE_Q;
   int kq_max = ranges[5] / SIZE_Q;
   for (int ky=ky_min; ky <= ky_max; ky++) {
    if (!ranges_overlap_bucket(&a->bucketY[ky], ranges)) continue;
    if (channel == 0) return 1;
    for (int ki=ki_min; ki <= ki_max; ki++) {
      if (!ranges_overlap_bucket(&a->bucketYI[ky][ki], ranges)) continue;
      if (channel == 1) return 1;
      for (int kq=kq_min; kq <= kq_max; kq++) {
        if (ranges_overlap_bucket(&a->bucketYIQ[ky][ki][kq],ranges)) return 1;
      }
    }
   }
   return 0;
}
bitvector compute_range(colors *a, int min, int max, int c, pixel_t *known, int qz, int guess) {
   if (c<0 || a->used==0 || a->max_diff[c] > 1) return make_vector(0);
   int size = max-min+1;
   bitvector v = make_vector(size);
   int ranges[6] = {0,0,0,255,510,510};
   for(int i=0; i<c; i++) {
     ranges[i] = UNSTRETCH(known->d[i]);
     ranges[3+i] = UNSTRETCH(known->d[i]);
   }
   for(int i=0; i<size; i++) {
     ranges[c] = clamp(UNSTRETCH((min+i) * qz + guess),0,(c==0?255:510));
     ranges[3+c] = clamp(UNSTRETCH((min+i) * qz + guess),0,(c==0?255:510));
     if (color_exists(a, &ranges[0],c)) set_bit(&v,i);
   }
   return v;
} */

static inline int pixel_in_bucket(const int y, const int i, const int q, const color_bucket *b, const unsigned int c) {
    int count = b->count;
    if (count > 0) {
         if (count < MAX_PER_BUCKET) {
             for (unsigned int j=0; j < count; j++) {
                 if (b->color[j].d[0] != y) continue;
                 if (b->color[j].d[1] != i) continue;
                 if (b->color[j].d[2] != q) continue;
                 return 1;
             }
         } else {
             if (b->color[1].d[0] < y) return 0;
             if (b->color[0].d[0] > y) return 0;
             if (b->color[1].d[1] < i) return 0;
             if (b->color[0].d[1] > i) return 0;
             if (b->color[1].d[2] < q) return 0;
             if (b->color[0].d[2] > q) return 0;
             return 1;
         }
    }
    return 0;
}

static inline int bucket2bv(const color_bucket *b, const int y, const int i, bitvector *bv) {
    int count = b->count;
    int nb = 0;
    if (count > 0) {
         if (count < MAX_PER_BUCKET) {
             for (unsigned int j=0; j < count; j++) {
                 if (b->color[j].d[0] != y) continue;
                 if (b->color[j].d[1] != i) continue;
                 bv->bits[ b->color[j].d[2] ] = 1;
                 nb++;
//                 fprintf(stderr,"exact match: %i\n",b->color[j].d[2]);
             }
         } else {
             if (b->color[1].d[0] < y) return 0;
             if (b->color[0].d[0] > y) return 0;
             if (b->color[1].d[1] < i) return 0;
             if (b->color[0].d[1] > i) return 0;
             for (int q = b->color[0].d[2]; q <=  b->color[1].d[2]; q++) {
                bv->bits[q] = 1;
                nb++;
//               fprintf(stderr,"range match: %i\n",q);
             }

         }
    }
    return nb;
}

static inline int c_r_2(colors *a, const unsigned int ty, const unsigned int ti, const int qz, const int guess, const int min, const int max, bitvector *bv) {
          unsigned int ky= ty/SIZE_Y;
          unsigned int ki= ti/SIZE_I;
          pixel_t p = {};
          p.d[0] = ty;
          p.d[1] = ti;
          unsigned int qmin = UNSTRETCH((min) * qz + guess);
          unsigned int qmax = UNSTRETCH((max) * qz + guess);
          unsigned int kqmin = qmin/SIZE_Q;
          unsigned int kqmax = qmax/SIZE_Q;
          int count = 0;
 
          for (int i=qmin; i<qmax; i++) {
                bv->bits[i]=0;
          }
//          fprintf(stderr,"counting...\n");
          for(unsigned int kq = kqmin; kq<=kqmax ; kq++) {
                count += bucket2bv(&a->bucketYIQ[ky][ki][kq],ty,ti,bv);
//                fprintf(stderr,"bucket %i,%i,%i -> count=%i\n",ky,ki,kq,count);
          }
          bv->pzero = &bv->bits[UNSTRETCH(guess)];
          bv->pzeroc = NULL;
//          fprintf(stderr,"counted %i\n",count);

/*          
          unsigned int tq = UNSTRETCH((min+0) * qz + guess);
          int count=0;
          for(unsigned int i=0; i<size; i++) {
//            p.d[2] = clamp(UNSTRETCH((min+i) * qz + guess),0,510);
//            int tq = UNSTRETCH((min+i) * qz + guess);
            assert(tq >= 0);
            assert(tq <= 510);
//            p.d[2] = tq;            
//            unsigned int kq = p.d[2]/SIZE_Q;
  //          if (pixel_in_bucket(ty,ti,tq,&a->bucketYIQ[ky][ki][tq/SIZE_Q],2)) set_bit(bv,i);
            int val = pixel_in_bucket(ty,ti,tq,&a->bucketYIQ[ky][ki][tq/SIZE_Q],2);
            change_bit(bv,i,val);
            count += val;
            tq++;
          }
          */
          
//          bv->pzero = &bv->bits[-min];
          return count;
}
static inline int c_r_1(colors *a, const unsigned int ty, const int qz, const int guess, const int min, const int max, bitvector *bv) {
//          unsigned int ti = UNSTRETCH((min+0) * qz + guess);

          unsigned int tizero = UNSTRETCH(guess);          
          bv->pzero = &a->bucketYI[ty][tizero];
          bv->pzeroc = &a->bucketYIc[ty][tizero];
          return count_range(bv,min,max);
/*          
          int count=0;
          for(unsigned int i=0; i<size; i++) {
//            unsigned int ti = clamp(UNSTRETCH((min+i) * qz + guess),0,510);
            //int ti = UNSTRETCH((min+i) * qz + guess);
            assert(ti >= 0);
            assert(ti <= 510);
//            if (a->bucketYI[ty][ti]) set_bit(bv,i);
            int val=a->bucketYI[ty][ti];
            change_bit(bv,i,val);
            count += val;
            ti++;
          }
          return count;
*/          
}
static inline int c_r_0(colors *a, const int qz, const int guess, const int min, const int max, bitvector *bv) {

//          unsigned int ty = UNSTRETCH((min+0) * qz + guess);
          unsigned int tyzero = UNSTRETCH(guess);
          bv->pzero = &a->bucketY[tyzero];
          bv->pzeroc = &a->bucketYc[tyzero];
          return count_range(bv,min,max);
//          return 2;
/*          
          int count=0;
          for(unsigned int i=0; i<size; i++) {
//            unsigned int ty = clamp(UNSTRETCH((min+i) * qz + guess),0,255);
            assert(ty >= 0);
            assert(ty <= 510);
//            if (ty < 0) continue;
//            if (ty > 255) break;
//            if (a->bucketY[ty]) set_bit(bv,i);
            int val = a->bucketY[ty];
            change_bit(bv,i,val);
            count += val;
            ty++; // = UNSTRETCH((min+i) * qz + guess);
          }
          return count;
          
          */
}

int compute_range(colors *a, const int min, const int max, const int c, const pixel_t *known, const int qz, const int guess, bitvector *bv) {
   if (c<0 || a->used==0 || a->max_diff[c] > 1) {init_vector(bv,0); return max-min+1;}
//   if (c<0 || a->used==0 ) {init_vector(0); return;}
   int size = max-min+1;
//   if (size < 2) {return NULL;}
//   init_vector(bv,size);
   init_vector_size(bv,size);
   bv->min = min;
   bv->max = max;

   switch (c) {
       case 0:
//          c_r_0(a,qz,guess,min,size, bv);
                // we're not doing this if quantization is > 1 so might as well optimize for that
          return c_r_0(a,COLOR_STRETCH,guess,min,max, bv);
          break;

       case 1:
//          c_r_1(a,UNSTRETCH(known->d[0]),qz,guess,min,size, bv);
          return c_r_1(a,UNSTRETCH(known->d[0]),COLOR_STRETCH,guess,min,max, bv);
          break;

       case 2:
//          c_r_2(a,UNSTRETCH(known->d[0]),UNSTRETCH(known->d[1]),qz,guess,min,size, bv);
          return c_r_2(a,UNSTRETCH(known->d[0]),UNSTRETCH(known->d[1]),COLOR_STRETCH,guess,min,max, bv);
          break;
   }
   return -1;
//   compute_cumul(bv,0,size-1);
}



void input_bounds(color_bucket *b, pixel_t *p, int n) {
  if (n==0)   for(int c=0; c<3; c++)   if (p->d[c] < b->color[0].d[c]) b->color[0].d[c] = p->d[c];
  if (n==1)   for(int c=0; c<3; c++)   if (p->d[c] > b->color[1].d[c]) b->color[1].d[c] = p->d[c];
}

void update_bounds(color_bucket *b, pixel_t *p) {
//  fprintf(stderr,"update_bounds: Y:%i..%i I:%i..%i Q:%i..%i    to ", b->color[0].d[0], b->color[1].d[0], b->color[0].d[1], b->color[1].d[1], b->color[0].d[2], b->color[1].d[2]);
  for(int c=0; c<3; c++) {
    if (p->d[c] < b->color[0].d[c]) b->color[0].d[c] = p->d[c];
    if (p->d[c] > b->color[1].d[c]) b->color[1].d[c] = p->d[c];
    assert(p->d[c] >= b->color[0].d[c]);
    assert(p->d[c] <= b->color[1].d[c]);
  }
//  fprintf(stderr,"Y:%i..%i I:%i..%i Q:%i..%i \n", b->color[0].d[0], b->color[1].d[0], b->color[0].d[1], b->color[1].d[1], b->color[0].d[2], b->color[1].d[2]);
}
void init_bounds(color_bucket *b) {
  pixel_t first = b->color[0];
  pixel_t second = b->color[1];
  update_bounds(b,&first);
  update_bounds(b,&second);
  for (int k=2; k < b->count; k++) {
    update_bounds(b,&b->color[k]);
  }
}


int auto_indexing_heuristic(colors *a) {
  int exact_colors = 0;
  int continuous_area = 0;
  int exact_buckets = 0;
  int continuous_buckets = 0;
  
  for (int ky=0; ky<BUCKETS_Y; ky++) {
   for (int ki=0; ki<BUCKETS_I; ki++) {
    for (int kq=0; kq<BUCKETS_Q; kq++) {
     int c = a->bucketYIQ[ky][ki][kq].count;
     if (c>0) {
        if (c < MAX_PER_BUCKET) {
                exact_buckets++;
                exact_colors += c;
        } else {
                continuous_buckets++;
                int yrange = a->bucketYIQ[ky][ki][kq].color[1].d[0] - a->bucketYIQ[ky][ki][kq].color[0].d[0] + 1;
                int irange = a->bucketYIQ[ky][ki][kq].color[1].d[1] - a->bucketYIQ[ky][ki][kq].color[0].d[1] + 1;
                int qrange = a->bucketYIQ[ky][ki][kq].color[1].d[2] - a->bucketYIQ[ky][ki][kq].color[0].d[2] + 1;
                continuous_area += yrange*irange*qrange;
        }
     }
    }
   }
  }
  fprintf(stderr,"Color buckets: %i exact (%i colors), %i continuous (area=%.2f%%).\n",exact_buckets,exact_colors,continuous_buckets,100.0*continuous_area/(256*511*511/2));

  if (continuous_area + exact_colors < 256*511*511/30) {fprintf(stderr,"Small total color area, so "); return 1;}
  if (continuous_area + exact_colors > 256*511*511/5) {fprintf(stderr,"Large total color area, so "); return 0;}

  if (exact_buckets+continuous_buckets < BUCKETS_Y*BUCKETS_I*BUCKETS_Q/10) {fprintf(stderr,"Few buckets, so "); return 1;}

  if (2*continuous_buckets + exact_colors < 1000) {fprintf(stderr,"Small color table, so "); return 1;}


  if (continuous_buckets < BUCKETS_Y*BUCKETS_I*BUCKETS_Q/10) {fprintf(stderr,"Few continuous buckets, so "); return 1;}
//  if (exact_colors < 1500) {fprintf(stderr,"Few exact colors, so "); return 1;}
  fprintf(stderr,"Large color table, "); return 0;

}




#define BUCKET_CONVERSION_MIN_COLORS 3

// if a bucket with indexed colors has a hull that covers less than X % of the bucket space, convert it to a full (continuous) bucket
#define BUCKET_CONVERSION_THRESHOLD 1000
#define BUCKET_CONVERSION_THRESHOLD_C 5
void convert_indexed_clusters_to_full_buckets(colors *a) {

//  return;
  
  int nb_converted = 0;
  int nb_colors = 0;
  for (int ky=0; ky<BUCKETS_Y; ky++) {
   for (int ki=0; ki<BUCKETS_I; ki++) {
    for (int kq=0; kq<BUCKETS_Q; kq++) {
     int c = a->bucketYIQ[ky][ki][kq].count;
     if (c>BUCKET_CONVERSION_MIN_COLORS && c<MAX_PER_BUCKET) {
        int min[3]={-1,-1,-1};
        int max[3]={-1,-1,-1};
        for (int i=0; i<c; i++) {
         for (int chan=0; chan<3; chan++) {
          int v = a->bucketYIQ[ky][ki][kq].color[i].d[chan];
          if (min[chan]==-1 || v<min[chan]) min[chan]=v;
          if (max[chan]==-1 || v>max[chan]) max[chan]=v;
         }
        }
        if ( (max[0]-min[0]+1)*(max[1]-min[1]+1)*(max[2]-min[2]+1)*1000 < SIZE_Y*SIZE_I*SIZE_Q * (BUCKET_CONVERSION_THRESHOLD + (c) * BUCKET_CONVERSION_THRESHOLD_C)) {
/*                fprintf(stderr,"Converting partial cluster (containing %i colors) to full cluster (%i,%i,%i) covering %g < %i promille.\n",c,ky,ki,kq,
                        1000.0*(max[0]-min[0]+1)*(max[1]-min[1]+1)*(max[2]-min[2]+1) / (SIZE_Y*SIZE_I*SIZE_Q), 
                        (BUCKET_CONVERSION_THRESHOLD + (c) * BUCKET_CONVERSION_THRESHOLD_C)
                        );*/
                init_bounds(&a->bucketYIQ[ky][ki][kq]);
                nb_converted++;
                nb_colors += c;
                a->bucketYIQ[ky][ki][kq].count = MAX_PER_BUCKET;
       }
     }
    }
   }
  }
  if (nb_converted>0)  fprintf(stderr,"Converted %i partial buckets (%i colors) to full buckets.\n",nb_converted,nb_colors);
}
    
    
void init_bounds_input(color_bucket *b, int c) {
  if (c == MAX_PER_BUCKET) return;
  pixel_t first = b->color[0];
  pixel_t second = b->color[1];
  if (c>0) update_bounds(b,&first);
  if (c>1) update_bounds(b,&second);
  for (int k=2; k < c; k++) {
    update_bounds(b,&b->color[k]);
  }
}

int same_color(pixel_t a, pixel_t b) {
   for (int c=0; c<3; c++)
       if (a.d[c] != b.d[c]) return 0;
   return 1;
}
int bigger_color(pixel_t a, pixel_t b) {
   for (int c=0; c<3; c++) {
       if (a.d[c] < b.d[c]) return 0;
       if (a.d[c] > b.d[c]) return 1;
   }
   assert(2<1);
   return 1;  // equal, should not happen
}

void add_color_to_bucket(color_bucket *b, pixel_t *p) {
  int c = b->count;
  if (c < MAX_PER_BUCKET) {
    int k=0;
    for (; k < c; k++) {
        if (same_color(b->color[k],*p)) return;
        if (bigger_color(b->color[k],*p)) break;
    }
    for (int i=c; i > k; i--) {
        b->color[i] = b->color[i-1];
    }
    b->color[k] = *p;
//        fprintf(stderr,"Bucket[%i][%i][%i]: has %i colors, adding color (%i,%i,%i)\n",ky,ki,kq,c,y,i,q);
    b->count = c+1;
    if (c+1 == MAX_PER_BUCKET)  init_bounds(b);
//    fprintf(stderr,"added color (%i,%i,%i)\n",p->d[0],p->d[1],p->d[2]);
  } else {
    update_bounds(b,p);
  }
}
void add_color(colors *a, pixel_t *p) {
  int y = UNSTRETCH(p->d[0]);
  int i = UNSTRETCH(p->d[1]);
  int q = UNSTRETCH(p->d[2]);
  pixel_t up = {};
  up.d[0] = y;
  up.d[1] = i;
  up.d[2] = q;
  int ky = y / SIZE_Y;
  int ki = i / SIZE_I;
  int kq = q / SIZE_Q; 
  if (! bucket_exists(ky,ki,kq) ) fprintf(stderr,"adding color (%i,%i,%i) in non-existing bucket (%i,%i,%i)\n",y,i,q,ky,ki,kq);
   add_color_to_bucket(&a->bucketYIQ[ky][ki][kq], &up);
//  add_color_to_bucket(&a->bucketYI[ky][ki], &up);
  a->bucketY[y]=1;
//  add_color_to_bucket(&a->bucket, &up);
//  if (!pixel_in_bucket(&up,&a->bucketYIQ[ky][ki][kq],2)) fprintf(stderr,"add_color problem: %i,%i,%i \n",y,i,q);
//  if (!reduce_range_Q(a, &q, &q, y, i)) fprintf(stderr,"add_color problem: %i,%i,%i \n",y,i,q);
}



void fill_bucketYI(colors *a, int glob_minq, int glob_maxq) {
  for (int ky=0; ky<BUCKETS_Y; ky++) {
   for (int ki=0; ki<BUCKETS_I; ki++) {
    for (int kq=0; kq<BUCKETS_Q; kq++) {
     int c = a->bucketYIQ[ky][ki][kq].count;
     if (c>0) {
      if (c<MAX_PER_BUCKET) {
        for (int i=0; i<c; i++) {
          int ty = a->bucketYIQ[ky][ki][kq].color[i].d[0];
          int ti = a->bucketYIQ[ky][ki][kq].color[i].d[1];
          a->bucketYI[ty][ti]=1;

//           fprintf(stderr,"setting bucketYI[%i][%i]=1\n",ty,ti);
        }
      } else {
        int ty_min = a->bucketYIQ[ky][ki][kq].color[0].d[0];
        int ti_min = a->bucketYIQ[ky][ki][kq].color[0].d[1];
        int tq_min = a->bucketYIQ[ky][ki][kq].color[0].d[2];
        int ty_max = a->bucketYIQ[ky][ki][kq].color[1].d[0];
        int ti_max = a->bucketYIQ[ky][ki][kq].color[1].d[1];
        int tq_max = a->bucketYIQ[ky][ki][kq].color[1].d[2];
        for (int ty=ty_min; ty<=ty_max; ty++) {
         for (int ti=ti_min; ti<=ti_max; ti++) {
            int minq;
            int maxq;
            if (ty<63) {
              minq=254-2*ty+(abs(ti-255)/2)*2;
              maxq=256+2*ty;
            } else if (ty>=192) {
              minq=255-2*(255-ty);
              maxq=255+2*(255-ty)-((1+abs(ti-255))/2)*2;
            } else {
              minq=MAX(1+(ty-128)*2,128-(ty-63)*2+(abs(ti-255)/2)*2);
              maxq=MIN(382+(ty-63)*2,383+(191-ty)*2-((1+abs(ti-255))/2)*2);
            }
            if (glob_minq >= minq) minq=glob_minq;
            if (glob_maxq <= maxq) maxq=glob_maxq;
            if (glob_minq >= maxq) maxq=glob_minq;
            if (glob_maxq <= minq) minq=glob_maxq;

            if (minq > tq_max) continue;
            if (maxq < tq_min) continue;
           
            a->bucketYI[ty][ti]=1;
//           fprintf(stderr,"setting bucketYI[%i][%i]=1\n",ty,ti);
         }
        }
      }
     }
    }
   }
  }
  // compute cumuls
  int ycount = 0;
  a->bucketYc[0] = 0;
  for (int y=0; y<256; y++) {
    ycount += a->bucketY[y];
    a->bucketYc[y+1] = ycount;
    int icount = 0;
    a->bucketYIc[y][0] = 0;
    for (int i=0; i<511; i++) {
      icount += a->bucketYI[y][i];
      a->bucketYIc[y][i+1] = icount;
    }
  }
  
}
