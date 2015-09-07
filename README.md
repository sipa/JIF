# JIF
JiF (or JPiF?) Image Format.

JiF is a lossless image format based on MANIAC compression.
MANIAC (Meta-Adaptive Near-zero Integer Arithmetic Coding) is a variant of
CABAC (context-adaptive binary arithmetic coding), where the contexts are nodes of decision trees
which are dynamically learned at encode time.

JiF outperforms PNG, FFV1, lossless WebP, and lossless JPEG2000 in terms of compression ratio.

Moreover, JiF supports a form of progressive interlacing (essentially a generalization/improvement of PNG's Adam7)
which means that any prefix (e.g. partial download) of a compressed file can be used as a reasonable
lossy encoding of the entire image.

JiF is now called FLIF and it moved over here: https://github.com/jonsneyers/FLIF
