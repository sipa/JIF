#include "yiq.h"

#include <stdlib.h>

ColorVal static inline get_min_y() {
    return 0;
}

ColorVal static inline get_max_y() {
    return 255;
}

ColorVal static inline get_min_i(ColorVal y) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());

    if (y<63) {
      return 252-4*y;
    } else if (y>=192) {
      return 3+4*(y-192);
    } else {
      return 0;
    }
}

ColorVal static inline get_max_i(ColorVal y) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());

    if (y<63) {
      return 258+4*y;
    } else if (y>=192) {
      return 507-4*(y-192);
    } else {
      return 510;
    }
}

ColorVal static inline get_min_q(ColorVal y, ColorVal i) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());
    assert(i >= get_min_i(y));
    assert(i <= get_max_i(y));

    if (y<63) {
      return 254-2*y+(abs(i-255)/2)*2;
    } else if (y>=192) {
      return 255-2*(255-y);
    } else {
      return std::max(1+(y-128)*2, 128-(y-63)*2+(abs(i-255)/2)*2);
    }
}

ColorVal static inline get_max_q(ColorVal y, ColorVal i) {
    assert(y >= get_min_y());
    assert(y <= get_max_y());
    assert(i >= get_min_i(y));
    assert(i <= get_max_i(y));

    if (y<63) {
      return 256+2*y;
    } else if (y>=192) {
      return 255+2*(255-y)-((1+abs(i-255))/2)*2;
    } else {
      return std::min(382+(y-63)*2, 383+(191-y)*2-((1+abs(i-255))/2)*2);
    }
}

ColorVal ColorRangesYIQ::min(int p, int r, int c) const {
    switch(p) {
    case 0:
        return get_min_y();
    case 1:
        return get_min_i((*image)(0,r,c));
    case 2:
        return get_min_q((*image)(0,r,c), (*image)(1,r,c));
    }
    assert(false);
    return -1;
}

ColorVal ColorRangesYIQ::max(int p, int r, int c) const {
    switch(p) {
    case 0:
        return get_max_y();
    case 1:
        return get_max_i((*image)(0,r,c));
    case 2:
        return get_max_q((*image)(0,r,c), (*image)(1,r,c));
    }
    assert(false);
    return -1;
}

