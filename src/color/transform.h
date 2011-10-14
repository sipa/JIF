//#include "plane.h"

// echte YIQ: y = 0.299*r + 0.587*g + 0.114*b
//            i = 0.595716*r - 0.274453*g -0.321263*b
//            q = 0.211456*r - 0.522591*g +0.311135*b

// benadering: y = 0.25*r + 0.5*g + 0.25*b
//             i = 1 * r          - 1 * b
//             q = 0.5*r  - 1*g   + 0.5*b

YIQA_Image RGBA_to_YIQA(RGBA_Image &in) {
        YIQA_Image out(in.height, in.width);
        out.A = in.A;
        for (int y=0; y<in.height; y++) {
          for (int x=0; x<in.width; x++) {
            int R = in.R.get(x,y);
            int G = in.G.get(x,y);
            int B = in.B.get(x,y);
            int Y = ((R + B) / 2 + G) / 2;
            int I = R - B;
            int Q = (R + B) / 2 - G;

            out.Y.set(x,y,Y);
            out.I.set(x,y,I);
            out.Q.set(x,y,Q);
          }
        }
        return out;
}

RGBA_Image YIQA_to_RGBA(YIQA_Image &in) {
       RGBA_Image out(in.height, in.width);
       out.A = in.A;
       for (int y=0; y<in.height; y++) {
          for (int x=0; x<in.width; x++) {
            int Y = in.Y.get(x,y);
            int I = in.I.get(x,y);
            int Q = in.Q.get(x,y);
            int R = Y + (Q + 257) / 2 + (I + 257) / 2 - 256;
            int G = Y - (Q + 256) / 2 + 128;
            int B = Y + (Q + 257) / 2 - (I + 256) / 2;

            out.R.set(x,y,R);
            out.G.set(x,y,G);
            out.B.set(x,y,B);
          }
       }
       return out;
}