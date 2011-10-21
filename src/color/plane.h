#include <vector>

#include <stdint.h>
#include <assert.h>

typedef uint16_t pixel_t;
static const pixel_t MIN_PIXEL_VAL = 0;
static const pixel_t MAX_PIXEL_VAL = 65535;

#define NB_PLANE_SEMANTICS 7
enum plane_semantics {RED, GREEN, BLUE, ALPHA, Y, I, Q};

// some common color spaces
static const plane_semantics aRGBA[4] = {RED, GREEN, BLUE, ALPHA};
static const std::vector<plane_semantics> RGBA(aRGBA, aRGBA+4);
static const plane_semantics aYIQA[4] = {Y, I, Q, ALPHA};
static const std::vector<plane_semantics> YIQA(aYIQA, aYIQA+4);


class Plane {
public:
  uint16_t height;      // layer size (could be smaller than canvas size)
  uint16_t width;

private:
  std::vector<std::vector<pixel_t> > p;

public:
  inline pixel_t get(uint16_t x, uint16_t y) {
        return p[y][x];
  }
  inline void set(uint16_t x, uint16_t y, pixel_t v) {
        assert(v >= MIN_PIXEL_VAL);
        assert(v <= MAX_PIXEL_VAL);
        p[y][x] = v;
  }
  Plane(uint16_t h, uint16_t w, pixel_t init_pixel_val = MIN_PIXEL_VAL) : height(h), width(w) {
        p.resize(h, std::vector<pixel_t>(w, init_pixel_val));
  }
};


class MetaData {
public:
  plane_semantics plane_type;
  uint8_t bit_depth;     // refers to the unstretched data
  uint16_t stretch;      // 1 = no stretch (stretch is used for more accurate interpolation)

  uint8_t zoom_factor_x; // 1 = no zoom
  uint8_t zoom_factor_y;

  uint32_t frame_number; // for animations/movies
  uint16_t z_index;      // determines order of overlapping layers (same frame, same plane type, lower index is "below" higher index)

  uint16_t offset_x;     //
  uint16_t offset_y;

  MetaData(plane_semantics pt, uint8_t bd = 8, uint16_t s = 1, uint8_t zx = 1, uint8_t zy = 1, uint32_t fn = 0, uint16_t zi = 0, uint16_t ox = 0, uint16_t oy = 0) :
        plane_type(pt), bit_depth(bd), stretch(s), zoom_factor_x(zx), zoom_factor_y(zy), frame_number(fn), z_index(zi), offset_x(ox), offset_y(oy) 
        {}
};



class ImageData {
public:
  std::vector<MetaData> info;
  std::vector<Plane> data;
  uint16_t height;  // canvas size
  uint16_t width;   // max res: 65536x65536  (4096 megapixels)

private:
  bool using_z_index;
  uint32_t nb_frames;
  std::vector<std::vector<int> > reverse_index;

  void construct_reverse_index() {
    reverse_index.clear();
    for (unsigned int p = 0; p < nb_frames; p++) {
        std::vector<int> ri(NB_PLANE_SEMANTICS, -1);
        for(unsigned int i=0; i < info.size(); i++) {
                if (info[i].frame_number == p)
                        ri[info[i].plane_type] = i;
        }
        reverse_index.push_back(ri);
    }
  }

public:
  Plane & plane(plane_semantics p, uint32_t frame = 0) {
        assert(frame <= nb_frames);
        assert( ! using_z_index); // todo: handle multiple (possibly overlapping) planes per frame
        assert(reverse_index[frame][p] >= 0); // -1 if the plane does not exist
        return data[reverse_index[frame][p]];
  }
  pixel_t inline get(plane_semantics p, uint16_t x, uint16_t y, uint32_t frame = 0) {
        assert(frame <= nb_frames);
        assert( ! using_z_index); // todo: handle multiple (possibly overlapping) planes per frame
        assert(reverse_index[frame][p] >= 0);
        return data[reverse_index[frame][p]].get(x,y);
  }
  void inline set(plane_semantics p, uint16_t x, uint16_t y, pixel_t pixel, uint32_t frame = 0) {
        assert(frame <= nb_frames);
        assert( ! using_z_index); // todo: handle multiple (possibly overlapping) planes per frame
        data[reverse_index[frame][p]].set(x,y,pixel);
  }
  void remove_plane(plane_semantics p, uint32_t frame = 0) {
        assert(frame <= nb_frames);
        assert( ! using_z_index); // todo: handle multiple (possibly overlapping) planes per frame
        assert(reverse_index[frame][p] >= 0);
        info.erase(info.begin()+reverse_index[frame][p]);
        data.erase(data.begin()+reverse_index[frame][p]);
        construct_reverse_index();
        assert(reverse_index[frame][p] == -1);
  }
  void add_plane(plane_semantics p, uint32_t frame = 0) {
        assert(frame <= nb_frames);
        assert( ! using_z_index); // todo: handle multiple (possibly overlapping) planes per frame
        assert(info.size() == data.size());
        int pos = info.size();
        info.push_back(MetaData(p));
        data.push_back(Plane(height, width));
        construct_reverse_index();
        assert(reverse_index[frame][p] == pos);
  }
  // constructor that takes image size and a list of plane types, and creates a one-frame image where each plane type fits the canvas size
  ImageData(uint16_t h, uint16_t w, std::vector<plane_semantics> pt = RGBA) : height(h), width(w) {
        using_z_index = false;
        nb_frames = 1;
        int nb_planes = pt.size();
        std::vector<int> ri(NB_PLANE_SEMANTICS, -1);
        for(int i = 0; i < nb_planes; i++) {
                info.push_back(MetaData(pt[i]));
                data.push_back(Plane(h, w));
                ri[pt[i]] = i;
        }
        reverse_index.push_back(ri);
  }
};


