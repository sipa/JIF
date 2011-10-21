//#include "plane.h"



// echte YIQ: y = 0.299*r + 0.587*g + 0.114*b
//            i = 0.595716*r - 0.274453*g -0.321263*b
//            q = 0.211456*r - 0.522591*g +0.311135*b

// benadering: y = 0.25*r + 0.5*g + 0.25*b
//             i = 1 * r          - 1 * b
//             q = 0.5*r  - 1*g   + 0.5*b

class RGBToBWTransformer : public Transformer {

  class RGBToBWTransformerData : public TransformerData { }

  void transform(const ImageData &in, ImageData &out) { 
    out = in;
    out.remove_plane(RED);
    out.remove_plane(GREEN);
    out.remove_plane(BLUE);
    Plane &RP = in.plane(RED);
    Plane &GP = in.plane(GREEN);
    Plane &BP = in.plane(BLUE);
    MetaData &RM = in.metadata(RED);
    MetaData &GM = in.metadata(GREEN);
    MetaData &BM = in.metadata(BLUE);
    if (RM.info->is_simple() && GM.info->is_simple() && BM.info->is_simple()) {
    } else {
      
    }

    out.add_plane(Y);
    Plane &YP = out.plane(Y);
    MetaData &YM = in.metadata(Y);
    
  }
}

ImageData RGB_to_YIQ(ImageData &in, uint32_t frame = 0) {
        ImageData out = in;
        out.remove_plane(RED, frame);
        out.remove_plane(GREEN, frame);
        out.remove_plane(BLUE, frame);
        out.add_plane(Y);
        out.add_plane(I);
        out.add_plane(Q);
        Plane &RP = in.plane(RED, frame);
        Plane &GP = in.plane(GREEN, frame);
        Plane &BP = in.plane(BLUE, frame);
        Plane &YP = out.plane(Y, frame);
        Plane &IP = out.plane(I, frame);
        Plane &QP = out.plane(Q, frame);
        for (int y=0; y<in.height; y++) {
          for (int x=0; x<in.width; x++) {
            int R = RP.get(x,y);
            int G = GP.get(x,y);
            int B = BP.get(x,y);

            int vY = ((R + B) / 2 + G) / 2;
            int vI = R - B + 255;
            int vQ = (R + B) / 2 - G + 255;

            YP.set(x,y,vY);
            IP.set(x,y,vI);
            QP.set(x,y,vQ);
          }
        }
        return out;
}

ImageData YIQ_to_RGB(ImageData &in, uint32_t frame = 0) {
        ImageData out = in;
        out.remove_plane(Y, frame);
        out.remove_plane(I, frame);
        out.remove_plane(Q, frame);
        out.add_plane(RED, frame);
        out.add_plane(GREEN, frame);
        out.add_plane(BLUE, frame);
        Plane &YP = in.plane(Y, frame);
        Plane &IP = in.plane(I, frame);
        Plane &QP = in.plane(Q, frame);
        Plane &RP = out.plane(RED, frame);
        Plane &GP = out.plane(GREEN, frame);
        Plane &BP = out.plane(BLUE, frame);
        for (int y=0; y<in.height; y++) {
          for (int x=0; x<in.width; x++) {
            int vY = YP.get(x,y);
            int vI = IP.get(x,y);
            int vQ = QP.get(x,y);

            int R = vY + (vQ + 2) / 2 + (vI + 2) / 2 - 256;
            int G = vY - (vQ + 1) / 2 + 128;
            int B = vY + (vQ + 2) / 2 - (vI + 1) / 2;

            RP.set(x,y,R);
            GP.set(x,y,G);
            BP.set(x,y,B);
          }
        }
        return out;
}