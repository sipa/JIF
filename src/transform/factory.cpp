#include <string>

#include "transform.h"

#include "yiq.h"
#include "bounds.h"

Transform *create_transform(std::string desc)
{
    if (desc == "YIQ")
        return new TransformYIQ();
    if (desc == "BND")
        return new TransformBounds();
    return NULL;
}
