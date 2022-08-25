/**
* @file         ImageCodegen.hpp
*
* Purpose       Define the Codegen context for the Image class
*
* Authors       Dinesh Iyer
*
* Copyright 2020 The MathWorks, Inc.
*
*/

#ifndef imagesdatatypesimagecodegen_header
#define imagesdatatypesimagecodegen_header

#ifndef IMAGE_IN_CODEGEN_MODE
#define IMAGE_IN_CODEGEN_MODE   // Macro that specifies that header is being called in generated code
#endif

// In codegen, this header files are in a "flat" directory structure.
#include "datatypes_util.hpp"
#include "ImageDefs.hpp"
#include "ImageMetadata.hpp"
#include "Image.hpp"
#include "ImageUtils.hpp"
#include "ImageExceptions.hpp"


#endif // imagesdatatypesimagecodegen_header
