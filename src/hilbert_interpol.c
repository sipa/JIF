#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "config.h"
#include "image/image.h"
#include "color.h"
#include "interpol.h"

int median(int a, int b, int c, int d) {
  int ab0 = (a>b?b:a);
//  int ab1 = (a>b?a:b);
  int cd0 = (c>d?d:c);
//  int cd1 = (c>d?c:d);
  int t1 = (ab0<cd0?cd0:ab0);
//  int t2 = (ab1<cd1?ab1:cd1);
  return t1;
//  return (t1+t2+1)/2;  
}
int median3(int a, int b, int c) {
  int ab0 = (a>b?b:a);
  int ab1 = (a>b?a:b);
  return (c>ab1?ab1:(c>ab0?c:ab0));
}

void static bilinear(pixel_t *out, image_t *img, int x, int y, int x1, int x2, int y1, int y2, int mask) {
  wpixel_t p={};
  //fprintf(stderr,"bil: (%i,%i) in (%i,%i)-(%i,%i)\n",x,y,x1,y1,x2,y2);
  wpixel_add(&p,image_pixel(img,x1,y1),(x2-x)*(y2-y),mask);
  wpixel_add(&p,image_pixel(img,x2,y1),(x-x1)*(y2-y),mask);
  wpixel_add(&p,image_pixel(img,x1,y2),(x2-x)*(y-y1),mask);
  wpixel_add(&p,image_pixel(img,x2,y2),(x-x1)*(y-y1),mask);
  wpixel_div(&p,(x2-x1)*(y2-y1),mask);
  wpixel_set(out,&p,mask);
}

void static linear(pixel_t *out, image_t *img, int x, int y, int x1, int x2, int y1, int y2, int mask) {
  pixel_linear(out,image_pixel(img,x1,y1),(x-x1)+(y-y1),image_pixel(img,x2,y2),(x2-x)+(y2-y),mask);
}


void static linear_w(int f, pixel_t *out, image_t *img, int x, int y, int x1, int x2, int y1, int y2, int mask) {
  pixel_linear_w(f,out,image_pixel(img,x1,y1),(x-x1)+(y-y1),image_pixel(img,x2,y2),(x2-x)+(y2-y),mask);
}


void interpolate_recurse(image_t *img, int x1, int x2, int y1, int y2, int mask, int hastop, int hasright, int hasbottom, int hasleft) {
//  return;
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
//  if (px1+1>=px2 && py1+1>=py2) return;
  if (x2-x1 == 0 && y2-y1 == 0) return;
  if (x2 - x1 > y2 - y1) {
    int xm = (x1 + x2) / 2;
    if (!hastop)    linear(image_pixel(img,xm,y1),img,xm,y1,px1,px2,py1,py1,mask);
    if (!hasbottom) linear(image_pixel(img,xm,y2),img,xm,y2,px1,px2,py2,py2,mask);
    interpolate_recurse(img,x1  ,xm,y1,y2,mask,hastop,0,hasbottom,hasleft);     // left
    interpolate_recurse(img,xm+1,x2,y1,y2,mask,hastop,hasright,hasbottom,1);    // right
  } else {
    int ym = (y1 + y2) / 2;
    if (!hasleft)  linear(image_pixel(img,x1,ym),img,x1,ym,px1,px1,py1,py2,mask);
    if (!hasright) linear(image_pixel(img,x2,ym),img,x2,ym,px2,px2,py1,py2,mask);
    interpolate_recurse(img,x1,x2,y1  ,ym,mask,hastop,hasright,0,hasleft);      // top
    interpolate_recurse(img,x1,x2,ym+1,y2,mask,1,hasright,hasbottom,hasleft);   // bottom
  }
  
}



void static quadratic(pixel_t *out, image_t *img, int x, int y, int x0, int x1, int x2, int y0, int y1, int y2, int mask) {
  if (x0 < 0) x0 = 0;
  if (y0 < 0) y0 = 0;
  if (x1 < 0) x1 = 0;
  if (y1 < 0) y1 = 0;
  if (x0==x1 && y0==y1) {
    pixel_linear(out,image_pixel(img,x1,y1),(x-x1)+(y-y1),image_pixel(img,x2,y2),(x2-x)+(y2-y),mask);
  } else {
    wpixel_t p={};
    wpixel_add(&p,image_pixel(img,x0,y0),2,mask);
    int d=2;
    if (x0 < x1) {
      wpixel_add(&p,image_pixel(img,x0+1,y0),1,mask);
      d++;
    }
    if (y0 < y1) {
      wpixel_add(&p,image_pixel(img,x0,y0+1),1,mask);
      d++;
    }
    if (x0 > 0) {
      wpixel_add(&p,image_pixel(img,x0-1,y0),1,mask);
      d++;
    } 
    if (y0 > 0) {
      wpixel_add(&p,image_pixel(img,x0,y0-1),1,mask);
      d++;
    }
    wpixel_div(&p,d,mask);
    pixel_t p0={.d= {0,0,0}};
    wpixel_set(&p0,&p,mask);
    
    pixel_quadratic(out,&p0,image_pixel(img,x1,y1),image_pixel(img,x2,y2),mask);
  }    
}
void static cubic(pixel_t *out, image_t *img, int x, int y, int x0, int x1, int x2, int x3, int y0, int y1, int y2, int y3, int mask) {
  if (x0 < 0) x0 = 0;
  if (y0 < 0) y0 = 0;
  if (x1 < 0) x1 = 0;
  if (y1 < 0) y1 = 0;
  if (x2 < 0) x2 = 0;
  if (y2 < 0) y2 = 0;
  if (x0==x1 && y0==y1) {
    pixel_linear(out,image_pixel(img,x2,y2),(x-x2)+(y-y2),image_pixel(img,x3,y3),(x3-x)+(y3-y),mask);
  } else {
    pixel_cubic(out,image_pixel(img,x0,y0),image_pixel(img,x1,y1),image_pixel(img,x2,y2),image_pixel(img,x3,y3),mask);
  }    
}

int static box_size(int ix1, int ix2, int iy1, int iy2) {
  return (ix2-ix1+iy2-iy1)*2;
}

void static box_getxyd(int ix1, int ix2, int iy1, int iy2, int pos, int *x, int *y, int *d) {
  while (pos<0) pos+=box_size(ix1,ix2,iy1,iy2);
  while (pos >= 0) {
    if (pos < ix2 - ix1) { // top
      *x = ix1 + pos;
      *y = iy1;
      *d = (pos > 0);
      return;
    }
    pos -= ix2 - ix1;
    if (pos < iy2 - iy1) { // right
      *x = ix2;
      *y = iy1 + pos;
      *d = 2 + (pos > 0);
      return;
    }
    pos -= iy2 - iy1;
    if (pos < ix2 - ix1) { // bottom
      *x = ix2 - pos;
      *y = iy2;
      *d = 4 + (pos > 0);
      return;
    }
    pos -= ix2 - ix1;
    if (pos < iy2 - iy1) { // left
      *x = ix1;
      *y = iy2 - pos;
      *d = 6 + (pos > 0);
      return;
    }
    pos -= iy2 - iy1;
  }
  assert(0);
  // prevent warning that output arguments may be used uninitialized
  *d=0;
  *x=0;
  *y=0;
}

static const int dx[]={ 0, 1, 0,-1, 0, 1, 0,-1};
static const int dy[]={-1, 0, 1, 0,-1, 0, 1, 0};

void interpolate_diffuse(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1=x1-(x1>0);
  int py1=y1-(y1>0);
  int px2=x2;
  int py2=y2;
  wpixel_t sum= {};
  int pcount=0;
  //uint8_t *cc=calloc((px2-px1+1)*(py2-py1+1),1);
  void upd(int x, int y, const wpixel_t *p) {
    assert(!((x==px1 && y==py1) || (x==px2 && y==py1) || (x==px1 && y==py2) || (x==px2 && y==py2)));
    assert(x>=x1 && x<=x2 && y>=y1 && y<=y2);
    //cc[(px2-px1+1)*(y-py1)+x-px1]=1;
    wpixel_set(image_pixel(img,x,y),p,mask);
    //fprintf(stderr,"upd (%i,%i) to [%i,%i,%i]\n",x,y,image_pixel(img,x,y)->d[0],image_pixel(img,x,y)->d[1],image_pixel(img,x,y)->d[2]);
  }
  void updl(int x, int y, const pixel_t *p1, int d1, const pixel_t *p2, int d2) {
    assert(!((x==px1 && y==py1) || (x==px2 && y==py1) || (x==px1 && y==py2) || (x==px2 && y==py2)));
    assert(x>=x1 && x<=x2 && y>=y1 && y<=y2);
    //cc[(px2-px1+1)*(y-py1)+x-px1]=1;
    pixel_linear(image_pixel(img,x,y),p1,d1,p2,d2,mask);
  }
  const pixel_t
    *tl=image_pixel(img,px1,py1),
    *tr=image_pixel(img,px2,py1),
    *bl=image_pixel(img,px1,py2),
    *br=image_pixel(img,px2,py2);
  // corner pixels of the ``P-box''
  wpixel_add(&sum,tl,1,mask); pcount++;
  wpixel_add(&sum,tr,1,mask); pcount++;
  wpixel_add(&sum,bl,1,mask); pcount++;
  wpixel_add(&sum,br,1,mask); pcount++;
  // lines on the ``P-box''
  for (int x=px1+1; x<px2; x++) {
    { // top
      if (y1==0) updl(x,py1,tl,x-px1,tr,px2-x);
      pixel_t *p=image_pixel(img,x,py1);
      wpixel_add(&sum,p,1,mask); pcount++;
    }
    { // bottom
      updl(x,py2,bl,x-px1,br,px2-x);
      pixel_t *p=image_pixel(img,x,py2);
      wpixel_add(&sum,p,1,mask); pcount++;
    }
  }
  for (int y=py1+1; y<py2; y++) {
    { // left
      if (x1==0) updl(px1,y,tl,y-py1,bl,py2-y);
      pixel_t *p=image_pixel(img,px1,y);
      wpixel_add(&sum,p,1,mask); pcount++;
    }
    { // right
      updl(px2,y,tr,y-py1,br,py2-y);
      pixel_t *p=image_pixel(img,px2,y);
      wpixel_add(&sum,p,1,mask); pcount++;
    }
  }
  // inner box
  int ix1=px1+1;
  int ix2=px2-1;
  int iy1=py1+1;
  int iy2=py2-1;
  wpixel_div(&sum,pcount,mask);
  /*for (int x=ix1; x<=ix2; x++) {
    for (int y=iy1; y<=iy2; y++) {
      pixel_set(image_pixel(coder->rec,x,y),&sum,mask);
    }
  }*/
  int levels=MIN((ix2-ix1+1)/2,(iy2-iy1+1)/2);
  while (ix1<ix2 && iy1<iy2) { // successively smaller bounding boxes
    int len=box_size(ix1,ix2,iy1,iy2);
    for (int i=0; i<len; i++) { // pixels on bounding box
      int x,y,d;
      box_getxyd(ix1,ix2,iy1,iy2,i,&x,&y,&d);
      wpixel_t pix={};
      int rd=d/2;
      int Ux=dx[rd],Dx=-Ux,Rx=dx[rd+1],Lx=-Rx;
      int Uy=dy[rd],Dy=-Uy,Ry=dy[rd+1],Ly=-Ry;
      if (!(d & 1)) {
        /* X is written as: (2*(p1+p3+p5+p8) + 3*(p2+p4+p6+p7))/20
         *
         *    p1 p2 p3
         *    p4 X  p6
         *    p5 p7 p8
         *
         *    - p1-p5 are known
         *    - p6-p8 are interpolated themselves
         */
        pixel_t *p1=image_pixel(img,x+Lx+Ux,y+Ly+Uy),
                *p2=image_pixel(img,x+Ux,y+Uy),
                *p3=image_pixel(img,x+Rx+Ux,y+Ry+Uy),
                *p4=image_pixel(img,x+Lx,y+Ly),
                *p5=image_pixel(img,x+Lx+Dx,y+Ly+Dy);
        pixel_t p6,p7,p8;
        bilinear(&p6,img,x+Rx   ,y+Ry   ,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        bilinear(&p7,img,x+Dx   ,y+Dy   ,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        bilinear(&p8,img,x+Rx+Dx,y+Ry+Dy,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        //bilinear(&p6,img,x+Rx   ,y+Ry   ,px1,px2,py1,py2,mask);
        //bilinear(&p7,img,x+Dx   ,y+Dy   ,px1,px2,py1,py2,mask);
        //bilinear(&p8,img,x+Rx+Dx,y+Ry+Dy,px1,px2,py1,py2,mask);
        wpixel_add(&pix, p1,3,mask);
        wpixel_add(&pix, p2,2,mask);
        wpixel_add(&pix, p3,3,mask);
        wpixel_add(&pix, p4,2,mask);
        wpixel_add(&pix, p5,3,mask);
        wpixel_add(&pix,&p6,2,mask);
        wpixel_add(&pix,&p7,2,mask);
        wpixel_add(&pix,&p8,3,mask);
        wpixel_div(&pix,20,mask);
      } else {
        /* X is written as (9*(p3+p8)+5*(p2+p7+p4+p9)+4*(p1+p6+p5+p10))/54
         *
         *   p1 p2 p3 p4 p5
         *         X
         *   p6 p7 p8 p9 p10
         */
        pixel_t *p1=image_pixel(img,x+Ux+2*Lx,y+Uy+2*Ly),
                *p2=image_pixel(img,x+Ux+  Lx,y+Uy+  Ly),
                *p3=image_pixel(img,x+Ux     ,y+Uy     ),
                *p4=image_pixel(img,x+Ux+  Rx,y+Uy+  Ry),
                *p5=image_pixel(img,x+Ux+2*Rx,y+Uy+2*Ry);
        pixel_t p6,p7,p8,p9,p10;
        bilinear(&p6 ,img,x+Dx+2*Lx,y+Dy+2*Ly,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        bilinear(&p7 ,img,x+Dx+  Lx,y+Dy+  Ly,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        bilinear(&p8 ,img,x+Dx     ,y+Dy     ,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        bilinear(&p9 ,img,x+Dx+  Rx,y+Dy+  Ry,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        bilinear(&p10,img,x+Dx+2*Rx,y+Dy+2*Ry,ix1-1,ix2+1,iy1-1,iy2+1,mask);
        wpixel_add(&pix,  p1,4,mask);
        wpixel_add(&pix,  p2,5,mask);
        wpixel_add(&pix,  p3,9,mask);
        wpixel_add(&pix,  p4,5,mask);
        wpixel_add(&pix,  p5,4,mask);
        wpixel_add(&pix, &p6,4*31,mask);
        wpixel_add(&pix, &p7,5*31,mask);
        wpixel_add(&pix, &p8,9*31,mask);
        wpixel_add(&pix, &p9,5*31,mask);
        wpixel_add(&pix,&p10,4*31,mask);
        wpixel_div(&pix,27*32,mask);
      }
      upd(x,y,&pix);
      //bilinear(image_pixel(img,x,y),img,x,y,px1,px2,py1,py2,mask);
    }
    ix1++;
    ix2--;
    iy1++;
    iy2--;
    levels--;
  }
  while (ix1<=ix2 && iy1<=iy2) {
    upd(ix1,iy1,&sum);
    if (ix1<ix2) ix1++; else if (iy1<iy2) iy1++; else break;
  }
  /*for (int x=x1; x<=x2; x++) {
    for (int y=y1; y<=y2; y++) {
      if (!((x==px1 && y==py1) || (x==px2 && y==py1) || (x==px1 && y==py2) || (x==px2 && y==py2))) {
        if (!cc[(px2-px1+1)*(y-py1)+x-px1]) fprintf(stderr,"(%i,%i) wasn't updated s=(%i..%i,%i..%i) p=(%i..%i,%i..%i)!\n",x,y,x1,x2,y1,y2,px1,px2,py1,py2);
      }
    }
  }
  free(cc); */
}

void interpolate_recurse_gradient(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  if (px1+1>=px2 && py1+1>=py2) return;
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
    if (y1==0) linear(image_pixel(img,xm,py1),img,xm,py1,px1,px2,py1,py1,mask);
    gradient(image_pixel(img,xm,py2),img,xm,py2,px1,px2,py1,py2,mask);
    interpolate_recurse_gradient(img,x1  ,xm,y1,y2,mask);
    interpolate_recurse_gradient(img,xm+1,x2,y1,y2,mask);
  } else {
    int ym = (py1 + py2) / 2;
    if (x1==0) linear(image_pixel(img,px1,ym),img,px1,ym,px1,px1,py1,py2,mask);
    gradient(image_pixel(img,px2,ym),img,px2,ym,px1,px2,py1,py2,mask);
    interpolate_recurse_gradient(img,x1,x2,y1  ,ym,mask);
    interpolate_recurse_gradient(img,x1,x2,ym+1,y2,mask);
  }
}


void interpolate_recurse_biquadratic(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  if (px1+1>=px2 && py1+1>=py2) return;
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
//    int ix1 = xm - (px2-xm);
    int ix1 = px1;
    int ix0 = ix1 - (px2-xm);
    if (y1==0) quadratic(image_pixel(img,xm,py1),img,xm,py1,ix0,ix1,px2,py1,py1,py1,mask);
               quadratic(image_pixel(img,xm,py2),img,xm,py2,ix0,ix1,px2,py2,py2,py2,mask);
    interpolate_recurse_biquadratic(img,x1  ,xm,y1,y2,mask);
    interpolate_recurse_biquadratic(img,xm+1,x2,y1,y2,mask);
  } else {
    int ym = (py1 + py2) / 2;
//    int iy1 = ym - (py2-ym);
    int iy1 = py1;
    int iy0 = iy1 - (py2-ym);
    if (x1==0) quadratic(image_pixel(img,px1,ym),img,px1,ym,px1,px1,px1,iy0,iy1,py2,mask);
               quadratic(image_pixel(img,px2,ym),img,px2,ym,px2,px2,px2,iy0,iy1,py2,mask);
    interpolate_recurse_biquadratic(img,x1,x2,y1  ,ym,mask);
    interpolate_recurse_biquadratic(img,x1,x2,ym+1,y2,mask);
  }
}

void interpolate_recurse_bicubic(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  if (px1+1>=px2 && py1+1>=py2) return;
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
//    int ix1 = xm - (px2-xm);
    int ix1 = px1;
    int ix0 = ix1 - (px2-xm);
    int ix0b = ix0 - (px2-xm);
    if (y1==0) cubic(image_pixel(img,xm,py1),img,xm,py1,ix0b,ix0,ix1,px2,py1,py1,py1,py1,mask);
               cubic(image_pixel(img,xm,py2),img,xm,py2,ix0b,ix0,ix1,px2,py2,py2,py2,py2,mask);
    interpolate_recurse_bicubic(img,x1  ,xm,y1,y2,mask);
    interpolate_recurse_bicubic(img,xm+1,x2,y1,y2,mask);
  } else {
    int ym = (py1 + py2) / 2;
//    int iy1 = ym - (py2-ym);
    int iy1 = py1;
    int iy0 = iy1 - (py2-ym);
    int iy0b = iy0 - (py2-ym);
    if (x1==0) cubic(image_pixel(img,px1,ym),img,px1,ym,px1,px1,px1,px1,iy0b,iy0,iy1,py2,mask);
               cubic(image_pixel(img,px2,ym),img,px2,ym,px2,px2,px2,px2,iy0b,iy0,iy1,py2,mask);
    interpolate_recurse_bicubic(img,x1,x2,y1  ,ym,mask);
    interpolate_recurse_bicubic(img,x1,x2,ym+1,y2,mask);
  }
}


static const int weights[] = {3,3,1,1,2};
void interpolate_recurse_clipart(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  if (px1+1>=px2 && py1+1>=py2) return;
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
    const pixel_t *inputs[5] = {image_pixel(img,px1,py2),image_pixel(img,px2,py2),image_pixel(img,px1,py1),image_pixel(img,px2,py1),image_pixel(img,xm,py1)};
    if (y1==0) {
            pixel_med(image_pixel(img,xm,py1), inputs, weights, 4, mask);
    }
    //bottom
    pixel_med(image_pixel(img,xm,py2), inputs, weights, 5, mask);
    
    interpolate_recurse_clipart(img,x1  ,xm,y1,y2,mask);
    interpolate_recurse_clipart(img,xm+1,x2,y1,y2,mask);
  } else {
    int ym = (py1 + py2) / 2;
    const pixel_t *inputs[5] = {image_pixel(img,px2,py1),image_pixel(img,px2,py2),image_pixel(img,px1,py1),image_pixel(img,px1,py2),image_pixel(img,px1,ym)};
    if (x1==0) {
            pixel_med(image_pixel(img,px1,ym), inputs, weights, 4, mask);
    }
    //right
    pixel_med(image_pixel(img,px2,ym), inputs, weights, 5, mask);
    interpolate_recurse_clipart(img,x1,x2,y1  ,ym,mask);
    interpolate_recurse_clipart(img,x1,x2,ym+1,y2,mask);
  }
}

static const int weights_avg[] = {5,5,1,1,2};
void interpolate_recurse_weighted_avg(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  if (px1+1>=px2 && py1+1>=py2) return;
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
    const pixel_t *inputs[5] = {image_pixel(img,px1,py2),image_pixel(img,px2,py2),image_pixel(img,px1,py1),image_pixel(img,px2,py1),image_pixel(img,xm,py1)};
    if (y1==0) {
            pixel_avg(image_pixel(img,xm,py1), inputs, weights_avg, 4, mask);
    }
    //bottom
    pixel_avg(image_pixel(img,xm,py2), inputs, weights_avg, 5, mask);
    
    interpolate_recurse_weighted_avg(img,x1  ,xm,y1,y2,mask);
    interpolate_recurse_weighted_avg(img,xm+1,x2,y1,y2,mask);
  } else {
    int ym = (py1 + py2) / 2;
    const pixel_t *inputs[5] = {image_pixel(img,px2,py1),image_pixel(img,px2,py2),image_pixel(img,px1,py1),image_pixel(img,px1,py2),image_pixel(img,px1,ym)};
    if (x1==0) {
            pixel_avg(image_pixel(img,px1,ym), inputs, weights_avg, 4, mask);
    }
    //right
    pixel_avg(image_pixel(img,px2,ym), inputs, weights_avg, 5, mask);
    interpolate_recurse_weighted_avg(img,x1,x2,y1  ,ym,mask);
    interpolate_recurse_weighted_avg(img,x1,x2,ym+1,y2,mask);
  }
}


void interpolate_recurse2(int xf, int yf, image_t *img, int x1, int x2, int y1, int y2, int mask) {
  // xf, yf : 0...16
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  if (px1+1>=px2 && py1+1>=py2) return;
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
    if (y1==0) linear_w(xf,image_pixel(img,xm,py1),img,xm,py1,px1,px2,py1,py1,mask);
    linear_w(xf,image_pixel(img,xm,py2),img,xm,py2,px1,px2,py2,py2,mask);
    interpolate_recurse2(xf,yf,img,x1  ,xm,y1,y2,mask);
    interpolate_recurse2(xf,yf,img,xm+1,x2,y1,y2,mask);
  } else {
    int ym = (py1 + py2) / 2;
    if (x1==0) linear_w(yf,image_pixel(img,px1,ym),img,px1,ym,px1,px1,py1,py2,mask);
    linear_w(yf,image_pixel(img,px2,ym),img,px2,ym,px2,px2,py1,py2,mask);
    interpolate_recurse2(xf,yf,img,x1,x2,y1  ,ym,mask);
    interpolate_recurse2(xf,yf,img,x1,x2,ym+1,y2,mask);
  }
}


/* interpolate a single pixel: if neither top pixels or left side pixels are known, use bilinear interpolation,
 otherwise using barycentric interpolation on the remaining corner points and projections on top/left side */
void static interpolate_pixel(image_t *img, int x, int y, int x1, int x2, int y1, int y2, pixel_t *out, int mask) {
  int px1=x1 - (x1>0),px2=x2;
  int py1=y1 - (y1>0),py2=y2;
  pixel_t *bc[3] = {};
  int bx[3] = {}, by[3] = {};
  switch ((x1>0) + 2 * (y1>0)) {
    case 0: { // bilinear interpolation
      uint64_t f11 = (px2 - x) * (py2 - y);
      uint64_t f21 = (x - px1) * (py2 - y);
      uint64_t f12 = (px2 - x) * (y - py1);
      uint64_t f22 = (x - px1) * (y - py1);
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          uint64_t q11 = image_pixel(img,px1,py1)->d[c] * f11;
          uint64_t q21 = image_pixel(img,px2,py1)->d[c] * f21;
          uint64_t q12 = image_pixel(img,px1,py2)->d[c] * f12;
          uint64_t q22 = image_pixel(img,px2,py2)->d[c] * f22;
          out->d[c] = (2 * (q11 + q21 + q12 + q22) + (px2 - px1) * (py2 - py1)) / ((px2 - px1) * (py2 - py1) * 2);
          if (//out->d[c] < 0 ||
          out->d[c] > (c == 0 ? 25500 : 51000)) {
            fprintf(stderr,"out->d[%i] = %i\n",c,out->d[c]);
          }
        }
      }
      return;
    }
    case 1: {
      bc[0] = image_pixel(img,px1,y);
      bx[0] = px1;
      by[0] = y;
      bc[1] = image_pixel(img,px2,py1);
      bx[1] = px2;
      by[1] = py1;
      bc[2] = image_pixel(img,px2,py2);
      bx[2] = px2;
      by[2] = py2;
      break;
    }
    case 2: {
      bc[0] = image_pixel(img,x,py1);
      bx[0] = x;
      by[0] = py1;
      bc[1] = image_pixel(img,px1,py2);
      bx[1] = px1;
      by[1] = py2;
      bc[2] = image_pixel(img,px2,py2);
      bx[2] = px2;
      by[2] = py2;
      break;
    }
    case 3: {
      bc[0] = image_pixel(img,px2,py2);
      bx[0] = px2;
      by[0] = py2;
      bc[1] = image_pixel(img,x,py1);
      bx[1] = x;
      by[1] = py1;
      bc[2] = image_pixel(img,px1,y);
      bx[2] = px1;
      by[2] = y;
      break;
    }
  }
  for (int n = 0; n < 3; n++) {
    if (bx[n] == x && by[n] == y) {
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) out->d[c] = bc[n]->d[c];
      }
      return;
    }
  }
  int64_t qd = ((by[2] - by[0]) * (bx[1] - bx[2])) + ((bx[0] - bx[2]) * (by[1] - by[2]));
  int64_t q0 = ((by[2] - y) * (bx[1] - bx[2])) + ((x - bx[2]) * (by[1] - by[2]));
  int64_t q1 = ((by[2] - by[0]) * (x - bx[2])) + ((bx[0] - bx[2]) * (y - by[2]));
  int64_t q2 = ((y - by[0]) * (bx[1] - x)) + ((bx[0] - x) * (by[1] - y));
  if (qd < 0) {
    qd = -qd;
    q0 = -q0;
    q1 = -q1;
    q2 = -q2;
  }
  //  printf("qd=%lli q0=%lli q1=%lli q2=%lli\n",(long long)qd,(long long)q0,(long long)q1,(long long)q2);
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = (2 * (q0 * bc[0]->d[c] + q1 * bc[1]->d[c] + q2 * bc[2]->d[c]) + qd - (qd > 3 ? 1 : 0)) / (2 * qd);
      if (//out->d[c] < 0 ||
      out->d[c] > (c == 0 ? 25500 : 51000)) {
        fprintf(stderr,"out->d[%i] = %i\n",c,out->d[c]);
      }
    }
  }
}

void interpolate_bilibary(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  for (int y = y1; y <= y2; y++) {
    for (int x = x1; x <= x2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
        interpolate_pixel(img,x,y,x1,x2,y1,y2,image_pixel(img,x,y),mask);
      }
    }
  }
}

void static interpolate_pixel_clipart(image_t *img, int x, int y, int x1, int x2, int y1, int y2, pixel_t *out, int mask) {
  int px1=x1 - (x1>0),px2=x2;
  int py1=y1 - (y1>0),py2=y2;
  int mx=(x2-x1)/2;
  int my=(y2-y1)/2;
  int m = (mx>my?my:mx);
  int p;
  
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          p = -1;
          if ( (x-x1 + y-y1) < m ) p = image_pixel(img,px1,py1)->d[c];
          if ( (x2-x + y-y1) < m ) p = image_pixel(img,px2,py1)->d[c];
          if ( (x-x1 + y2-y) < m ) p = image_pixel(img,px1,py2)->d[c];
          if ( (x2-x + y2-y) < m ) p = image_pixel(img,px2,py2)->d[c];
          if (p == -1) {
                int q11 = image_pixel(img,px1,py1)->d[c];
                int q21 = image_pixel(img,px2,py1)->d[c];
                int q12 = image_pixel(img,px1,py2)->d[c];
                int q22 = image_pixel(img,px2,py2)->d[c];
                p = median(q11,q21,q12,q22);
          }
          out->d[c] = p;
          if (out->d[c] > (c == 0 ? 25500 : 51000)) {
            fprintf(stderr,"out->d[%i] = %i\n",c,out->d[c]);
          }
        }
      }
      return;
}


void interpolate_clipart(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  for (int y = y1; y <= y2; y++) {
    for (int x = x1; x <= x2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
        interpolate_pixel_clipart(img,x,y,x1,x2,y1,y2,image_pixel(img,x,y),mask);
      }
    }
  }
}
void static interpolate_pixel_clipart2(image_t *img, int x, int y, int x1, int x2, int y1, int y2, pixel_t *out, int mask) {
  int px1=x1 - (x1>0),px2=x2;
  int py1=y1 - (y1>0),py2=y2;
  int xm=(x2+x1)/2;
  int ym=(y2+y1)/2;
  int mx=(x2-x1)/2;
  int my=(y2-y1)/2;
  int m = (mx>my?my:mx);
  pixel_t *top, *left, *bottomright;
  bottomright = image_pixel(img,px2,py2);
  
  if (x>xm) { top = image_pixel(img,px2,py1); }
       else { top = image_pixel(img,px1,py1); }
  if (y>ym) { left = image_pixel(img,px1,py2); }
       else { left = image_pixel(img,px1,py1); }

  switch ((x1>0) + 2 * (y1>0)) {
    case 0: { 
      break;
    }
    case 1: {
      left = image_pixel(img,px1,y);
      break;
    }
    case 2: {
      top = image_pixel(img,x,py1);
      break;
    }
    case 3: {
      top = image_pixel(img,x,py1);
      left = image_pixel(img,px1,y);
      break;
    }
  }
  pixel_t *pixel=bottomright;
  if ( (x-x1) < m && y-y1-(x-x1) > 0 ) pixel = left;
  if ( (y-y1) < m && y-y1-(x-x1) <= 0) pixel = top;
  
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          out->d[c] = pixel->d[c];
        }
      }
      return;
}


void interpolate_clipart2(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  for (int y = y1; y <= y2; y++) {
    for (int x = x1; x <= x2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
        interpolate_pixel_clipart2(img,x,y,x1,x2,y1,y2,image_pixel(img,x,y),mask);
      }
    }
  }
}
double dist(int x1, int y1, int x2, int y2) {
        return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}
void static interpolate_pixel_clipart3(image_t *img, int x, int y, int x1, int x2, int y1, int y2, pixel_t *out, int mask) {
  int px1=x1 - (x1>0),px2=x2;
  int py1=y1 - (y1>0),py2=y2;
  int xm=(x2+x1)/2;
  int ym=(y2+y1)/2;
  pixel_t *top, *left, *bottomright, *topleft;
  bottomright = image_pixel(img,px2,py2);

  if (y>0) { top = image_pixel(img,x,y-1); } 
        else {
          if (x>xm) { top = image_pixel(img,px2,py1); }
               else { top = image_pixel(img,px1,py1); }
        }
  if (x>0) { left = image_pixel(img,x-1,y); }
        else {
          if (y>ym) { left = image_pixel(img,px1,py2); }
               else { left = image_pixel(img,px1,py1); }
        }

  assert(x>0 || y>0);

  if (x==0) {
        topleft = left;
  } else {
        if (y==0) {
                topleft = top;
        } else {
                topleft = image_pixel(img,x-1,y-1);
        }
  }
        
  int include_bottomright = 0;
  
  if ( dist(x,y,x2,y2) <= dist(x1,y1,x2,y2)/2 ) include_bottomright = 1;
  if ( dist(x,y,x2,y2) <= dist(x1,y1,x2,y2)/4 ) include_bottomright = 2;
  
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          if (include_bottomright == 2) {
                out->d[c] = bottomright->d[c];
          }
          else if (include_bottomright == 1) {
                out->d[c] = median(top->d[c],left->d[c],topleft->d[c],bottomright->d[c]);
          } else {
                out->d[c] = median3(top->d[c],left->d[c],topleft->d[c]);
          }
        }
      }
      return;
}
void interpolate_clipart3(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  for (int y = y1; y <= y2; y++) {
    for (int x = x1; x <= x2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
        interpolate_pixel_clipart3(img,x,y,x1,x2,y1,y2,image_pixel(img,x,y),mask);
      }
    }
  }
}

void interpolate_nn(image_t *img, int x1, int x2, int y1, int y2, int mask) {
  int px1 = x1 - (x1>0);
  int py1 = y1 - (y1>0);
  int px2 = x2;
  int py2 = y2;
  int x,y;
  for (y = y1; y <= (y1+y2)/2; y++) {
    for (x = x1; x <= (x1+x2)/2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          image_pixel(img,x,y)->d[c] = image_pixel(img,px1,py1)->d[c];
        }
      }
      }
    }
    for ( ; x <= x2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          image_pixel(img,x,y)->d[c] = image_pixel(img,px2,py1)->d[c];
        }
      }
      }
    }
  }
  for ( ; y <= y2; y++) {
    for (x = x1; x <= (x1+x2)/2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          image_pixel(img,x,y)->d[c] = image_pixel(img,px1,py2)->d[c];
        }
      }
      }
    }
    for ( ; x <= x2; x++) {
      if (!((x == px1 && y==py1) || (x == px2 && y==py1) || (x == px1 && y==py2) || (x==px2 && y==py2))) {
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c))) {
          image_pixel(img,x,y)->d[c] = image_pixel(img,px2,py2)->d[c];
        }
      }
      }
    }
  }
}

static const char *interpolation_names_array[] = {
  "BiliBary",
  "Diffuse",
  "Recursive Bilinear",
  "Clipart (corners)",
  "Clipart (nearest known)",
  "Clipart (nearest corner)",
  "Clipart (median_of_3/bottomright circle)",
  "Recursive weighted median",
  "Recursive weighted average",
  "Recursive Biquadratic",
  "Recursive Bicubic",
  "Recursive Gradient",
  "Irrelevant (all corners identical)",
  NULL
};

const char **interpolation_name = interpolation_names_array;

void interpolation(int type, image_t *img, int x1, int x2, int y1, int y2, int mask, int hastop, int hasright, int hasbottom, int hasleft) {
  switch (type) {
        case 0: {
                interpolate_bilibary(img, x1, x2, y1, y2, mask);
                return;
                }
        case 1: {
                interpolate_diffuse(img, x1, x2, y1, y2, mask);
                return;
                }
        case 2: {
                interpolate_recurse(img, x1, x2, y1, y2, mask, hastop, hasright, hasbottom, hasleft);
                return;
                }
        case 3: {
                interpolate_clipart(img, x1, x2, y1, y2, mask);
                return;
                }
        case 4: {
                interpolate_clipart2(img, x1, x2, y1, y2, mask);
                return;
                }
        case 5: {
                interpolate_nn(img, x1, x2, y1, y2, mask);
                return;
                }
        case 6: {
                interpolate_clipart3(img, x1, x2, y1, y2, mask);
                return;
                }
        case 7: {
                interpolate_recurse_clipart(img, x1, x2, y1, y2, mask);
                return;
                }
        case 8: {
                interpolate_recurse_weighted_avg(img, x1, x2, y1, y2, mask);
                return;
                }
        case 9: {
                interpolate_recurse_biquadratic(img, x1, x2, y1, y2, mask);
                return;
                }
        case 10: {
                interpolate_recurse_bicubic(img, x1, x2, y1, y2, mask);
                return;
                }
        case 11: {
                interpolate_recurse_gradient(img, x1, x2, y1, y2, mask);
                return;
                }

/*        case 6: {
                interpolate_recurse2(4,4,img, x1, x2, y1, y2, mask);
                return;
        }
        case 7: {
                interpolate_recurse2(4,8,img, x1, x2, y1, y2, mask);
                return;
        }
        case 8: {
                interpolate_recurse2(4,12,img, x1, x2, y1, y2, mask);
                return;
        }
        case 9: {
                interpolate_recurse2(8,4,img, x1, x2, y1, y2, mask);
                return;
        }
        case 10: {
                interpolate_recurse2(12,8,img, x1, x2, y1, y2, mask);
                return;
        }
        case 11: {
                interpolate_recurse2(12,4,img, x1, x2, y1, y2, mask);
                return;
        }
        case 12: {
                interpolate_recurse2(12,12,img, x1, x2, y1, y2, mask);
                return;
        }
        */
  }
}
/*
void interpolate(image_t *img, int x1, int x2, int y1, int y2, int mask) {
        interpolation(INTERPOLATION_METHOD, img, x1, x2, y1, y2, mask);
}
*/
