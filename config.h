#ifndef _CONFIG_H_
#define _CONFIG_H_ 1

#define SYMBOL_STATS 0
#define USE_AVERAGE_CHANCES 0
#define MULTISCALE_CHANCES 1
#define MULTISCALE_LEVELS 4

#define INTERPOLATION_METHOD 2
#define ACCURATE_FRAC 0
#define ACCURATE_YIQ_RANGES 1

#define MIN_LONG_WIDTH 4
#define MIN_LONG_PIXELS 128
#define MIN_ASSYM_WIDTH 16
#define MIN_ASSYM_PIXELS 512

#define COLOR_STRETCH 128
#define DEPTH_LEVELS 32

#define USE_BORDER_GUESSING 0

#define LINE_FACTOR 1
#define MAX_AVG_ERROR 3
#define MAX_AVG_ERROR_LINE 2.5

/*
// context quantization
#define QZ_SIZE 16
#define QZ_LRDIFF 8
#define QZ_GRAD 1
#define QZ_CROSS 1
#define QZ_VAR 1
#define QZ_BVAR 10
#define QZ_SIZE_S 62
#define QZ_BVAR_S 10
*/

// context quantization
#define QZ_SIZE 16
#define QZ_SIZE_S 60
#define QZ_SIZE_SL 35

#define QZ_DEPTH 1
#define QZ_DEPTH_S 1
#define QZ_DEPTH_SL 1

#define QZ_LRDIFF 4
#define QZ_LRDIFF_S 6
#define QZ_GRAD 1
#define QZ_CROSS 1
#define QZ_VAR 1
#define QZ_BVAR 12
#define QZ_BVAR_S 5


#define QZ_Y 12
#define QZ_I 14

// not bigger than QZ_Y
#define QZ_YQ 7
#define QZ_YY 7

// not bigger than QZ_I
#define QZ_II 4

#define RAC_BITS 40

#endif
