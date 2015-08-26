#include <string>

#include "transform.h"

#include "yiq.h"
#include "bounds.h"
#include "zoom.h"
#include "plane_permute.h"
#include "subgreen.h"
#include "quantize.h"
#include "bitshuffle.h"
#include "colorbuckets.h"

Transform *create_transform(std::string desc)
{
    if (desc == "YIQ")
        return new TransformYIQ();
    if (desc == "ZOOM")
        return new TransformZoom();
    if (desc == "BND")
        return new TransformBounds();
    if (desc == "PLP")
        return new TransformPP();
    if (desc == "SGR")
        return new TransformSG();
    if (desc == "QTZ")
        return new TransformQuantize();
    if (desc == "BSH")
        return new TransformBS();
    if (desc == "ACB")
        return new TransformCB();
    return NULL;
}
