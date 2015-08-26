#ifndef _ZOOM_H_
#define _ZOOM_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"

class TransformZoom : public Transform {
#ifdef SMOOTHZOOM
protected:
    int zoomlevel;
public:
    void configure(const int setting) { zoomlevel = setting; }
    void data(Image& image) const {
      printf("Transform:Zoom\n");
      for (int z = 0; z <= image.zooms(); z++) {
      if (z % 2 == 0) {
//        fprintf(stdout,"Zoomlevel %i: horizontal scan, current size: %i rows, %i cols \t",z,image.rows(z),image.cols(z));
        // horizontal: scan the odd rows
        for (int p = 0; p < image.numPlanes(); p++) {
          for (int r = 1; r < image.rows(z); r += 2) {
            for (int c = 0; c < image.cols(z); c++) {
                    ColorVal top = image(p,z,r-1,c);
                    ColorVal bottom = image(p,z,r,c);
                    ColorVal avg = (top+bottom)/2;
                    image(p,z,r-1,c) = avg; //top+bottom;
                    image(p,z,r,c) = bottom + PARITYBIT*(1&(top+bottom));
//                    if (p==0) printf("image(0,%i,%i -1,%i) = top: %i, bottom: %i, avg: %i, rest:%i\n",z,r,c,top,bottom,avg, image(p,z,r,c));
            }
          }
        }
      } else {
//        fprintf(stdout,"Zoomlevel %i: vertical scan, current size: %i rows, %i cols \t",z,image.rows(z),image.cols(z));
        // vertical: scan the odd columns
        for (int p = 0; p < image.numPlanes(); p++) {
          for (int r = 0; r < image.rows(z); r++) {
            for (int c = 1; c < image.cols(z); c += 2) {
                    ColorVal left = image(p,z,r,c-1);
                    ColorVal right = image(p,z,r,c);
                    ColorVal avg = (left+right)/2;
                    image(p,z,r,c-1) = avg; //left+right;
                    image(p,z,r,c) = right + PARITYBIT*(1&(left+right));
//                    if (p==0) printf("image(0,%i,%i,%i -1) = left: %i, right: %i, avg: %i, rest:%i\n",z,r,c,left,right,avg, image(p,z,r,c));
            }
          }
        }
      }
      }
    }

    void invData(Image& image) const {
      printf("Transform:UnZoom\n");
      for (int z = image.zooms(); z >= zoomlevel; z--) {
//      for (int z = 0; z <= image.zooms(); z++) {
      if (z % 2 == 0) {
        fprintf(stdout,"Zoomlevel %i: horizontal scan, current size: %i rows, %i cols \n",z,image.rows(z),image.cols(z));
        // horizontal: scan the odd rows
        for (int p = 0; p < image.numPlanes(); p++) {
          for (int r = 1; r < image.rows(z); r += 2) {
            for (int c = 0; c < image.cols(z); c++) {
                    ColorVal avg = image(p,z,r-1,c);
                    ColorVal bottom = image(p,z,r,c);
                    ColorVal second = (bottom & PARITYMASK);
                    ColorVal top = 2*avg - second + ((bottom&PARITYBIT) ? 1 : 0);
                    image(p,z,r-1,c) = top;
                    image(p,z,r,c) = second;
//                    printf ("plane:%i, zoomlevel:%i, r:%i, c:%i --> ",p,z,r,c);
//                    printf("top=%i, avg=%i, bottom=%i, parity=%i\t",top,avg,second,(bottom&PARITYBIT));
//                    if (top < 0 || top > 255) printf("OOPS!");
//                    printf("\n");
            }
          }
        }
      } else {
        fprintf(stdout,"Zoomlevel %i: vertical scan, current size: %i rows, %i cols \n",z,image.rows(z),image.cols(z));
        // vertical: scan the odd columns
        for (int p = 0; p < image.numPlanes(); p++) {
          for (int r = 0; r < image.rows(z); r++) {
            for (int c = 1; c < image.cols(z); c += 2) {
                    ColorVal avg = image(p,z,r,c-1);
                    ColorVal right = image(p,z,r,c);
                    ColorVal second = (right & PARITYMASK);
                    ColorVal left = 2*avg - second + ((right&PARITYBIT) ? 1 : 0);
                    image(p,z,r,c-1) = left;
                    image(p,z,r,c) = second;
//                    printf ("plane:%i, zoomlevel:%i, r:%i, c:%i --> ",p,z,r,c);
//                    printf("left=%i, avg=%i, right=%i, parity=%i\t",left,avg,second,(right&PARITYBIT));
//                    if (left < 0 || left > 255) printf("OOPS!");
//                    printf("\n");
            }
          }
        }
      }
      }


      for (int z = zoomlevel-1; z >= 0; z--) {
      if (z % 2 == 0) {
        fprintf(stdout,"INTERPOLATING Zoomlevel %i: horizontal scan, current size: %i rows, %i cols \n",z,image.rows(z),image.cols(z));
        // horizontal: scan the odd rows
        for (int p = 0; p < image.numPlanes(); p++) {
          for (int r = 1; r < image.rows(z); r += 2) {
            for (int c = 0; c < image.cols(z); c++) {
                    ColorVal avg = image(p,z,r-1,c);
                    //if (r+1 < image.rows(z)) image(p,z,r,c) = (avg + image(p,z,r+1,c))/2;
                    //else 
                    image(p,z,r,c) = avg;
            }
          }
        }
      } else {
        fprintf(stdout,"INTERPOLATING Zoomlevel %i: vertical scan, current size: %i rows, %i cols \n",z,image.rows(z),image.cols(z));
        // vertical: scan the odd columns
        for (int p = 0; p < image.numPlanes(); p++) {
          for (int r = 0; r < image.rows(z); r++) {
            for (int c = 1; c < image.cols(z); c += 2) {
                    ColorVal avg = image(p,z,r,c-1);
                    //if (c+1 < image.cols(z)) image(p,z,r,c) = (avg + image(p,z,r,c+1))/2;
                    //else 
                    image(p,z,r,c) = avg;
            }
          }
        }
      }
      }

   }
#endif
};

#endif
