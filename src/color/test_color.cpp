#include <string>
#include <stdio.h>
#include <assert.h>

#include "plane.h"
#include "transform.h"


int main() {
  ImageData testRGB(1024,768);

  for (int y = 0; y < testRGB.height ; y++)
    for (int x = 0; x < testRGB.width ; x++) {
      int check = (x+y) % 256;
      testRGB.set(RED,x,y,check);
      testRGB.set(GREEN,x,y,x % 256);
      testRGB.set(BLUE,x,y,y % 256);
      testRGB.set(ALPHA,x,y,x % 256);
      if (testRGB.get(RED,x,y) != check) fprintf(stdout,"Problem: %i != %i\n",check,testRGB.get(RED,x,y));
    }


  ImageData testYIQ = RGB_to_YIQ(testRGB);

  ImageData testRGB2 = YIQ_to_RGB(testYIQ);

  bool same_image = true;
  int pixels_ok = 0;
  for (int y = 0; y < testRGB.height ; y++)
    for (int x = 0; x < testRGB.width ; x++) {
      if (testRGB.get(RED,x,y) != testRGB2.get(RED,x,y)) {same_image = false;
           fprintf(stdout,"x=%i, y=%i, red1=%i, red2=%i\n", x, y, testRGB.get(RED,x,y), testRGB2.get(RED,x,y));}
      if (testRGB.get(GREEN,x,y) != testRGB2.get(GREEN,x,y)) {same_image = false;
           fprintf(stdout,"x=%i, y=%i, g1=%i, g2=%i\n", x, y, testRGB.get(GREEN,x,y), testRGB2.get(GREEN,x,y));}
      if (testRGB.get(BLUE,x,y) != testRGB2.get(BLUE,x,y)) {same_image = false;
           fprintf(stdout,"x=%i, y=%i, b1=%i, b2=%i\n", x, y, testRGB.get(BLUE,x,y), testRGB2.get(BLUE,x,y));}
      if (testRGB.get(ALPHA,x,y) != testRGB2.get(ALPHA,x,y)) {same_image = false;
           fprintf(stdout,"x=%i, y=%i, a1=%i, a2=%i\n", x, y, testRGB.get(ALPHA,x,y), testRGB2.get(ALPHA,x,y));}

      pixels_ok++;
    }


  if (same_image) {
        fprintf(stdout,"All %i pixels OK",pixels_ok);
  } else {
        fprintf(stdout,"%i pixels OK, the others are not OK!", pixels_ok);
  }

  return 0;
}



