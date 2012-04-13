#!/bin/sh

find ../src -type f -a \( -name '*.cpp' -o -name '*.c' -o -name '*.h' \) -print0 | xargs -0 astyle --style=k/r -c --mode=c "$@"


