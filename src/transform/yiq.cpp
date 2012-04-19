#include "yiq.h"

#include <stdlib.h>

ColorVal static inline get_min_y(int par) {
    return 0;
}

ColorVal static inline get_max_y(int par) {
    return par*4-1;
}

ColorVal static inline get_min_i(int par, ColorVal y) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());

    if (y<par-1) {
      return 4*par-4-4*y;
    } else if (y>=3*par) {
      return 3+4*(y-3*par);
    } else {
      return 0;
    }
}

ColorVal static inline get_max_i(int par, ColorVal y) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());

    if (y<par-1) {
      return 4*par+2+4*y;
    } else if (y>=3*par) {
      return 8*par-5-4*(y-3*par);
    } else {
      return 8*par-2;
    }
}

ColorVal static inline get_min_q(int par, ColorVal y, ColorVal i) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());
    assert(i >= get_min_i(y));
    assert(i <= get_max_i(y));

    if (y<par-1) {
      return 4*par-2-2*y+(abs(i-4*par+1)/2)*2;
    } else if (y>=3*par) {
      return 4*par-1-2*(4*par-1-y);
    } else {
      return std::max(1+(y-2*par)*2, 2*par-(y-par+1)*2+(abs(i-4*par+1)/2)*2);
    }
}

ColorVal static inline get_max_q(int par, ColorVal y, ColorVal i) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());
    assert(i >= get_min_i(y));
    assert(i <= get_max_i(y));

    if (y<par-1) {
      return 4*par+2*y;
    } else if (y>=3*par) {
      return 4*par-1+2*(4*par-1-y)-((1+abs(i-4*par+1))/2)*2;
    } else {
      return std::min(6*par-2+(y-par+1)*2, 6*par-1+(3*par-1-y)*2-((1+abs(i-4*par+1))/2)*2);
    }
}

ColorVal ColorRangesYIQ::min(int p, int r, int c) const {
    switch(p) {
    case 0:
        return get_min_y(par);
    case 1:
        return get_min_i(par, (*image)(0,r,c));
    case 2:
        return get_min_q(par, (*image)(0,r,c), (*image)(1,r,c));
    }
    assert(false);
    return -1;
}

ColorVal ColorRangesYIQ::max(int p, int r, int c) const {
    switch(p) {
    case 0:
        return get_max_y(par);
    case 1:
        return get_max_i(par, (*image)(0,r,c));
    case 2:
        return get_max_q(par, (*image)(0,r,c), (*image)(1,r,c));
    }
    assert(false);
    return -1;
}

