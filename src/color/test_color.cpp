#include <string>
#include <stdio.h>
#include <assert.h>

#include "plane.h"
#include "transform.h"


int main() {
  RGBA_Image testRGB(1024,768);

  for (int y = 0; y < testRGB.height ; y++)
    for (int x = 0; x < testRGB.width ; x++) {
      testRGB.R.set(x,y,x+y % 256);
      testRGB.G.set(x,y,x % 256);
      testRGB.B.set(x,y,y % 256);
      testRGB.A.set(x,y,x % 256);
    }


  YIQA_Image testYIQ = RGBA_to_YIQA(testRGB);

  RGBA_Image testRGB2 = YIQA_to_RGBA(testYIQ);

  bool same_image = true;
  int pixels_ok = 0;
  for (int y = 0; y < testRGB.height ; y++)
    for (int x = 0; x < testRGB.width ; x++) {
      if (testRGB.R.get(x,y) != testRGB2.R.get(x,y)) same_image = false;
      if (testRGB.G.get(x,y) != testRGB2.G.get(x,y)) same_image = false;
      if (testRGB.B.get(x,y) != testRGB2.B.get(x,y)) same_image = false;
      if (testRGB.A.get(x,y) != testRGB2.A.get(x,y)) same_image = false;
      pixels_ok++;
    }


  if (same_image) {
        fprintf(stdout,"All %i pixels OK",pixels_ok);
  } else {
        fprintf(stdout,"Not OK!");
  }

  return 0;
}



