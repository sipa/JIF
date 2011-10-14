#include <vector>

#include <stdint.h>
#include <assert.h>

class P_8bit {
public:
  typedef uint8_t pixel_t;
  static const pixel_t MIN_PIXEL_VAL = 0;
  static const pixel_t MAX_PIXEL_VAL = 255;
};

class P_9bit {
public:
  typedef int16_t pixel_t;
  static const pixel_t MIN_PIXEL_VAL = -255;
  static const pixel_t MAX_PIXEL_VAL = 255;
};

template <typename Config> class Plane {
public:
  typedef typename Config::pixel_t pixel_t;

private:
  uint16_t height;      // max res: 65536x65536  (4096 megapixels)
  uint16_t width;
  std::vector<std::vector<pixel_t> > p;

public:
  pixel_t get(uint16_t x, uint16_t y) {
        return p[y][x];
  }
  void set(uint16_t x, uint16_t y, pixel_t v) {
        assert(v >= Config::MIN_PIXEL_VAL);
        assert(v <= Config::MAX_PIXEL_VAL);
        p[y][x] = v;
  }
  Plane(uint16_t h, uint16_t w, pixel_t init_pixel_val = Config::MIN_PIXEL_VAL) : height(h), width(w) {
        p.resize(h, std::vector<pixel_t>(w, init_pixel_val));
  }
};

class RGBA_Image {
public:
        Plane<P_8bit> R;
        Plane<P_8bit> G;
        Plane<P_8bit> B;
        Plane<P_8bit> A;
        uint16_t height;
        uint16_t width;

  RGBA_Image(uint16_t h, uint16_t w) : R(h,w), G(h,w), B(h,w), A(h,w), height(h), width(w)
  { }
};

class YIQA_Image {
public:
        Plane<P_8bit> Y;
        Plane<P_9bit> I;
        Plane<P_9bit> Q;
        Plane<P_8bit> A;
        uint16_t height;
        uint16_t width;

  YIQA_Image(uint16_t h, uint16_t w) : Y(h,w), I(h,w), Q(h,w), A(h,w), height(h), width(w)
  { }
};
