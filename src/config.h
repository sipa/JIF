#ifndef _CONFIG_H_
#define _CONFIG_H_ 1

#define SYMBOL_STATS 0
#define USE_AVERAGE_CHANCES 0
#define MULTISCALE_CHANCES 1
#define MULTISCALE_LEVELS 6

#define INTERPOLATION_METHOD 2
#define ACCURATE_FRAC 0
#define ACCURATE_YIQ_RANGES 1

#define MIN_LONG_WIDTH 4
#define MIN_LONG_PIXELS 128
#define MIN_ASSYM_WIDTH 16
#define MIN_ASSYM_PIXELS 512

#define COLOR_STRETCH 128
#define DEPTH_LEVELS 32
#define SPLIT_MODULO 3

// add extra quantization between 0 and 20/this number (roughly)
#define EXTRA_QUANTIZATION_DEPTH 4
#define EXTRA_QUANTIZATION_VAR 3

#define MAX_QUANTIZATION 30


#define PROGRESS_BAR_VERBOSITY 6
#define PROGRESS_BAR_CHAR "="

#define STRETCH(x) ((x)*COLOR_STRETCH+(COLOR_STRETCH/2))
#define UNSTRETCH(x) (((x)/COLOR_STRETCH))
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif


#define USE_BORDER_GUESSING 0
#define ALWAYS_SPLIT_HEURISTIC 0


#define SPLIT_CHANNELS_SEPARATELY 1

#define HILBERT_CURVE 0
#define FULL_BORDERS 0
#define FFV1 1


#define LINE_FACTOR 1
//#define MAX_AVG_ERROR 3
#define MAX_AVG_ERROR 50
#define MAX_AVG_ERROR_LINE 20


#define MAX_X 65535
#define MAX_Y 65535

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

#define QZ_LOG (32*10+9)

// context quantization
#define QZ_SIZE_S QZ_LOG
#define QZ_SIZE_SL QZ_LOG
#define QZ_DEPTH_S 1
#define QZ_DEPTH_SL 1
#define QZ_LRDIFF_S 1000
#define QZ_BVAR_S QZ_LOG


#define QZ_DEPTH 1

//#define QZ_SIZE 200
#define QZ_SIZE QZ_LOG
#define QZ_LRDIFF 1021
#define SIGNED_LRDIFF 1
#define QZ_GRAD QZ_LOG
#define QZ_CROSS QZ_LOG
#define QZ_VAR QZ_LOG
#define QZ_BVAR QZ_LOG
#define QZ_BDIFFE QZ_LOG
#define QZ_YG 256
#define QZ_IG 511
#define QZ_QG 511
#define QZ_YD 511
#define QZ_ID 1021



#define RAC_BITS 16

#define MAX_INTERPOLATIONS 12
#define MAX_INTERPOLATIONS_EVER 500


#define MAX_CONTEXTS   20000
#define MAX_CONTEXTS_1 5000
#define MAX_CONTEXTS_S 500
#define MAX_CONTEXT_SIZE_DEPTH_FACTOR 1
#define CONTEXT_VIRTUAL_OUTPUT_DELAY 4
#define CONTEXT_SPLIT_MIN_SIZE 20
#define CONTEXT_SPLIT_MAX_SIZE 200000
#define CONTEXT_SPLIT_MIN_SIZEP 20
#define CONTEXT_SPLIT_MAX_SIZEP 200000

#define CONTEXT_SPLIT_MIN_SIZE_S 15
#define CONTEXT_SPLIT_IMPROVE_FACTOR 1
//#define CONTEXT_SPLIT_IMPROVE_CONSTANT 5461*10
//#define CONTEXT_SPLIT_IMPROVE_CONSTANT_S 5461*8
#define CONTEXT_SPLIT_IMPROVE_CONSTANT 5461*10
#define CONTEXT_SPLIT_IMPROVE_CONSTANT_S 5461*8
                                // n bit improvement needed before splitting
#define CONTEXT_USE_IMPROVE_FACTOR 1
//#define CONTEXT_USE_IMPROVE_CONSTANT 5461*4
//#define CONTEXT_USE_IMPROVE_CONSTANT_S 5461*2
#define CONTEXT_USE_IMPROVE_CONSTANT 5461*3
#define CONTEXT_USE_IMPROVE_CONSTANT_S 5461*2

#define IGNORE_SPLIT_OPTIONS 0
#define CONTEXT_SPLIT_IGNORE_FACTOR 3
#define CONTEXT_SPLIT_IGNORE_COUNT 50

// if endpoint distance is bigger than this, always split (unit=quantization)
#define ALWAYS_SPLIT_LINE_DIST 3
// if border variance is bigger than this, always split (unit=quantization)
#define ALWAYS_SPLIT_AREA_DIST 1000000

#endif
