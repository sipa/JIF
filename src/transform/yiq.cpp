#include "yiq.h"

#include <stdlib.h>

/*
#ifndef SMOOTHZOOM
#ifndef PERMUTEPLANES

// this only works for pixels that actually correspond to real pixels (not averaged pixels)

ColorVal ColorRangesYIQ::min(int p, int r, int c) const {
    switch(p) {
    case 0:
        return get_min_y(par);
    case 1:
        return get_min_i(par, (*image)(0,r,c));
    case 2:
        return get_min_q(par, (*image)(0,r,c), (*image)(1,r,c));
    default:
        return min(p);
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
    default:
        return max(p);
    }
    assert(false);
    return -1;
}

ColorVal ColorRangesYIQ::min(int p, int z, int r, int c) const {
    switch(p) {
    case 0:
        return get_min_y(par);
    case 1:
        return get_min_i(par, (*image)(0,z,r,c));
    case 2:
        return get_min_q(par, (*image)(0,z,r,c), (*image)(1,z,r,c));
    default:
        return min(p);
    }
    assert(false);
    return -1;
}

ColorVal ColorRangesYIQ::max(int p, int z, int r, int c) const {
    switch(p) {
    case 0:
        return get_max_y(par);
    case 1:
        return get_max_i(par, (*image)(0,z,r,c));
    case 2:
        return get_max_q(par, (*image)(0,z,r,c), (*image)(1,z,r,c));
    default:
        return max(p);
    }
    assert(false);
    return -1;
}
#endif
#endif
*/