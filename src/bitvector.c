
#include "bitvector.h"


//static bitvector v = {};
//static int saved_nb = 0;


/*
void disable_vector() {
  saved_nb = v.nb;
  v.nb = -1;
}
void enable_vector() {
  v.nb = saved_nb;
}
*/

int inline get_bit(const bitvector *v, const int pos) {
  assert(pos <= v->nb);
  assert(pos >= 0);
  return v->bits[pos];
//  return v->cumul[pos+1] - v->cumul[pos];
}




int inline get_first_in_range(const bitvector *v, int min, int max) {
   if (min==max) return min;
   if (EMPTY_VECTOR(v)) return min;
   if (max > v->max) max = v->max;
   if (min < v->min) min = v->min;
   for (int i=min ; i <= max ; i++) {
        if (v->pzero[i]) { return i;}
   }
   
/*   
   int begin=min-orig_min;
   int end=max-orig_min;
   if (begin<0) begin=0;
   if (end>v->nb) end=v->nb;
   for (unsigned int i=begin ; i <= end ; i++) {
      if (get_bit(v,i)) return i+orig_min;
   }
   */
   assert(0);
   return -1;
}




