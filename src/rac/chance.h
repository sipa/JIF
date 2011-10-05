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
  BitChance zero;
  BitChance sign;
  BitChance exponent[MAX_BITS_PER_SYMBOL-1];
  BitChance mant[MAX_BITS_PER_SYMBOL-1];
}

class SymbolWriter {
  RacOutput* output;
  SymbolChance chance;

  SymbolWriter operator<<(int data) {

  void write_int(SymbolChance c, int(*range_test)(int, int), int *value, const int min, const int max) {

    assert(min<=max);
    assert(*value>=min);
    assert(*value<=max);
    assert(range_test(min,max));

    // avoid doing anything if the value is already known (this line is optional; behavior should be identical if commented out)
    if (range_test(min,max) == 1) return;

    if (*value) { // value is nonzero
      // only output zero bit if value could also have been zero
      if (max >= 0 && min <= 0 && range_test(0,0)) output.write(c.zero, false);
      bool sign = (*value > 0 ? true : false);
      // only output sign bit if value can be both pos and neg
      if (max > 0 && min < 0 && range_test(min,-1) && range_test(1,max)) output.write(c.sign, sign);
      if (sign && min <= 0) min = 1;
      if (!sign && max >= 0) max = -1;
      const unsigned int a = abs(*value);
      const unsigned int e = log2_tab[a];
      unsigned int emin = log2_tab[(sign ? abs(min) : abs(max))];
      unsigned int emax = log2_tab[(sign ? abs(max) : abs(min))];
      unsigned int i = emin;
      while (i < emax) {
        // if exponent >i is impossible, we are done
        if (sign && !range_test(1<<(i+1),max)) break;
        if (!sign && !range_test(min,-(1<<(i+1)))) break;
        // if exponent i is possible, output the exponent bit
        if ( sign && range_test(1<<i,1<<(i+1)-1)
         || !sign && !range_test(1-(1<<(i+1)),-(1<<i)))
                output.write(c.exponent[i], i == e);
        if (i == e) break;
        i++;
      }
      int have = (1 << e);
      int left = (1 << e)-1;
      for (unsigned int pos = e; pos>0;) {
        int bit = 1;
        int minval = (sign ? have : -(have+left));
        int maxval = (sign ? have+left : -have);
        if (min > minval) minval = min;
        if (max < maxval) maxval = max;
        left ^= (1 << (--pos));
        if (range_test(minval,maxval)==1) return;
        int minval1 = (sign ? have+(1<<pos) : minval);
        int maxval1 = (sign ? maxval : -(have+(1<<pos)));
        int minval0 = (sign ? minval : -(have+left));
        int maxval0 = (sign ? have+left : maxval);
        if (!range_test(minval1,maxval1)) { // 1-bit is impossible
           bit = 0;
        } else if (range_test(minval0,maxval0)) { // 0-bit and 1-bit are both possible
           bit = (a >> pos) & 1;
           output.write(c.mantissa[pos],bit);
        }
        have |= (bit << pos);
      }

    } else { // value is zero
      // only output zero bit if value could also have been nonzero
      if (range_test(min,-1) || range_test(1,max)) output.write(c.zero, true);
    }
  }
}


/* example:

  FILE* file = fopen("data.dat","w");
  RacOutput rac(file);
  
  SymbolWriter pixDiff[3] = { SymbolWriter(rac), SymbolWriter(rac), SymbolWriter(rac) };

  pixDiff[0] << pix1R << pix2R << pix3R;
  pixDiff[1] << pix1G << pix2G << pix3G;

*/
