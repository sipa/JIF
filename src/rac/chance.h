class BitChance {
  uint16_t chance;
}

class BitWriter {
  RacOutput* output;
  BitChance chance;
  
  BitWriter operator<<(bool data) {
  }
}

class SymbolChance {
  BitChance sign[9];
  BitChance expo[8];
  BitChance mant[9];
  BitChance zero;
}

class SymbolWriter {
  RacOutput* output;
  SymbolChance chance;

  SymbolWriter operator<<(int data) {
  }
}


/* example:

  FILE* file = fopen("data.dat","w");
  RacOutput rac(file);
  
  SymbolWriter pixDiff[3] = { SymbolWriter(rac), SymbolWriter(rac), SymbolWriter(rac) };

  pixDiff[0] << pix1R << pix2R << pix3R;
  pixDiff[1] << pix1G << pix2G << pix3G;

*/
