#ifndef _INTERPOL_H_
#define _INTERPOL_H_

#include "image/image.h"

#define MAX_INTERPOLATIONS 12
#define MAX_INTERPOLATIONS_EVER 500

extern const char **interpolation_name;

void interpolation(int type, image_t *img, int x1, int x2, int y1, int y2, int mask, int hastop, int hasright, int hasbottom, int hasleft);
//void interpolate(image_t *img, int x1, int x2, int y1, int y2, int mask);


void static gradient(pixel_t *out, image_t *img, int x, int y, int x1, int x2, int y1, int y2, int mask) {
   pixel_t *top, *left;

   top=image_pixel(img,x ,y1);
   if (y1==0) { pixel_linear(top,image_pixel(img,x1,y1),x - x1,image_pixel(img,x2,y1),x2 - x,mask); }
   if (y==y1) {out = top; return;}

   left=image_pixel(img,x1 ,y);
   if (x1==0) { pixel_linear(left,image_pixel(img,x1,y1),y - y1,image_pixel(img,x1,y2),y2 - y,mask); }
   if (x==x1) {out = left; return;}

   if (y==y2) {
        wpixel_t p={};
        wpixel_add(&p,image_pixel(img,x1,y1),-1,mask);
        wpixel_add(&p,image_pixel(img,x2,y1),-1,mask);
        wpixel_add(&p,image_pixel(img,x1,y2), 1,mask);
        wpixel_add(&p,image_pixel(img,x2,y2), 1,mask);
        wpixel_add(&p,top, 2,mask);
        wpixel_div(&p,2,mask);
        wpixel_clamp(&p,image_pixel(img,x1,y2),image_pixel(img,x2,y2),mask);
        wpixel_set(out,&p,mask);
        return;
  } else if (x==x2) {
        wpixel_t p={};
        wpixel_add(&p,image_pixel(img,x1,y1),-1,mask);
        wpixel_add(&p,image_pixel(img,x2,y1), 1,mask);
        wpixel_add(&p,image_pixel(img,x1,y2),-1,mask);
        wpixel_add(&p,image_pixel(img,x2,y2), 1,mask);
        wpixel_add(&p,left, 2,mask);
        wpixel_div(&p,2,mask);
        wpixel_clamp(&p,image_pixel(img,x2,y1),image_pixel(img,x2,y2),mask);
        wpixel_set(out,&p,mask);        
        return;
  }
  fprintf(stderr,"Should not happen.\n");

}

#endif

