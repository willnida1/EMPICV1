/**
* @file         datatypes_util.h
*
* Purpose       export/import switch for images::datatypes::Image classes
*
* Authors       Dinesh Iyer
*
* Copyright     2019-2020 The MathWorks, Inc.
*
*/
#ifdef SUPPORTS_PRAGMA_ONCE
#pragma once
#endif

#ifndef _IMAGES_DATATYPE_IMAGE_HPP_
#define _IMAGES_DATATYPE_IMAGE_HPP_

#ifdef IMAGE_IN_CODEGEN_MODE // Specifies this header is being used in generate code
#include "rtwtypes.h"

#ifndef IMAGESDATATYPESIMAGE_API
#define IMAGESDATATYPESIMAGE_API
#endif

#else // When not in generated code

#ifdef BUILDING_IMAGESDATATYPESIMAGE
/* For DLL_EXPORT_SYM and EXTERN_C */
#include "package.h"
/* For size_t, boolean_T, etc */
#include "tmwtypes.h"

#define IMAGESDATATYPESIMAGE_API  DLL_EXPORT_SYM

#else
// #define IMAGESDATATYPESIMAGE_API  DLL_IMPORT_SYM


#ifndef EXTERN_C
#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C extern
#endif
#endif

#ifndef IMAGESDATATYPESIMAGE_API
#define IMAGESDATATYPESIMAGE_API
#endif

#endif
#endif

#endif // _IMAGES_DATATYPE_IMAGE_HPP_
