#ifndef _INTERPOL_H_
#define _INTERPOL_H_

#include "image/image.h"
#include "indexing.h"


extern const char **interpolation_name;

void interpolation(int type, image_t *img, int x1, int x2, int y1, int y2, int mask, int hastop, int hasright, int hasbottom, int hasleft, colors *colors);
//void interpolate(image_t *img, int x1, int x2, int y1, int y2, int mask);

inline static void inv_dist_onepixel(pixel_t *pixel, image_t *img,int x,int y,int px1,int px2,int py1,int py2,int x1,int x2,int y1,int y2,int wx1,int wx2,int wy1,int wy2,int mask) {

        // do not touch the corners
        if (x == px1 && y == py1) return;
        if (x == px1 && y == py2) return;
        if (x == px2 && y == py1) return;
        if (x == px2 && y == py2) return;
        
        xwpixel_t p={};
        int64_t tweight = 0;
//        fprintf(stderr,"Pixel (%i,%i): ",y,x);
        for (int iy = wy1; iy <= wy2 ; iy++) {
          for (int ix = wx1; ix <= wx2 ; ix++) {
            if (iy >= y1 && iy <= y2 && ix >= x1 && ix <= x2
                    && (ix != px1 || iy != py1)
                    && (ix != px1 || iy != py2)
                    && (ix != px2 || iy != py1)
                    && (ix != px2 || iy != py2) ) continue;  // only consider known points
            if (ix == x && iy == y) continue; 
            int64_t weight = 1000000000/((x-ix)*(x-ix)+(y-iy)*(y-iy));
//            int weight = 100/(abs(x-ix)+abs(y-iy));
//            if (ix == px2 && iy == py2) weight *= 5;  // more weight for known corner
//            if (ix >= x2) weight /= 2;  // less weight for roughly approximated points
//            if (iy >= y2) weight /= 2;  // less weight for roughly approximated points
            
            xwpixel_addu(&p,image_pixel(img,ix,iy),weight,mask);
            tweight += weight;
//            fprintf(stderr,"building pixel (%i,%i) guess: using (%i,%i), total Y=%lli, I=%lli, Q=%lli\n",x,y,ix,iy,p.wd[0],p.wd[1],p.wd[2]);
            if ((tweight >> 55) > 0) fprintf(stderr,"overflow danger\n");
            if (tweight < 0) fprintf(stderr,"OVERFLOW\n");
          }
        }
//        fprintf(stderr,"tweight: %i \n",tweight);
        xwpixel_divs(&p,tweight,mask);
        xwpixel_set(pixel,&p,mask);
//        fprintf(stderr,"pixel (%i,%i) guess: Y=%i, I=%i, Q=%i\n",x,y,UNSTRETCH(pixel->d[0]),UNSTRETCH(pixel->d[1]),UNSTRETCH(pixel->d[2]));
}



void gradient(pixel_t *out, image_t *img, int x, int y, int x1, int x2, int y1, int y2, int mask);

#endif

