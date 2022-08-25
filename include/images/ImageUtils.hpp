/**
* @file         ImageUtils.hpp
*
* Purpose       Utility functions 
*
* Authors       Dinesh Iyer
*
* Copyright     2019-2020 MathWorks, Inc.
*
*/

#ifndef imagesdatatypesimage_utils_header
#define imagesdatatypesimage_utils_header

#ifdef IMAGE_IN_CODEGEN_MODE
#include "ImageDefs.hpp" /* To deal with the fact that PackNGo has no include file hierarchy during test */
#else
#include "images/datatypes/ImageDefs.hpp"
#endif

#include <cstring>
#include <string>
#include <map>

template<typename MapType, typename InputType>
static typename MapType::mapped_type findValueInMap(const MapType& map, const InputType& key) {
    typename MapType::const_iterator it = map.find(key);
    return it->second;
}

namespace images 
{
    namespace datatypes 
    {
        namespace utils 
        {
            static images::datatypes::ClassUnderlying getClassUnderlyingFromName(const std::string& className) 
            {
                using namespace images::datatypes;

                typedef std::map<std::string, ClassUnderlying> ClassMap;
                
                static const ClassMap classMap {
                    {"logical", ClassUnderlying::Logical},
                    {"uint8", ClassUnderlying::Uint8},
                    {"int8", ClassUnderlying::Int8},
                    {"uint16", ClassUnderlying::Uint16},
                    {"int16", ClassUnderlying::Int16},
                    {"uint32", ClassUnderlying::Uint32},
                    {"int32", ClassUnderlying::Int32},
                    {"single", ClassUnderlying::Single},
                    {"double", ClassUnderlying::Double} };

                    return findValueInMap(classMap, className);
            }

            static std::string getClassUnderlyingName(const images::datatypes::ClassUnderlying dtype) 
            {

                using namespace images::datatypes;

                typedef std::map<ClassUnderlying, std::string> DatatypeNameMap;

                static const DatatypeNameMap datatypeNameMap {
                    {ClassUnderlying::Logical, "logical"},
                    {ClassUnderlying::Uint8, "uint8"},
                    {ClassUnderlying::Int8, "int8"},
                    {ClassUnderlying::Uint16, "uint16"},
                    {ClassUnderlying::Int16, "int16"},
                    {ClassUnderlying::Uint32, "uint32"},
                    {ClassUnderlying::Int32, "int32"},
                    {ClassUnderlying::Single, "single"},
                    {ClassUnderlying::Double, "double"} };

                    return findValueInMap(datatypeNameMap, dtype);
            }

            static uint8_t getClassUnderlyingValue(const images::datatypes::ClassUnderlying dtype)
            {
                return static_cast<uint8_t>( dtype );
            }

            static images::datatypes::ColorFormat getColorFormatFromName(const std::string& pixelLayout) 
            {
                using namespace images::datatypes;

                typedef std::map<std::string, ColorFormat> ColorFormatMap;
                
                static const ColorFormatMap formatMap {
                    {"BGRA", ColorFormat::BGRA},
                    {"RGB", ColorFormat::RGB},
                    {"RGBA", ColorFormat::RGBA},
                    {"Grayscale",  ColorFormat::Grayscale},
                    {"CMYK",  ColorFormat::CMYK},
                    {"Spectral", ColorFormat::Spectral},
                    {"Lab", ColorFormat::Lab},
                    {"HSV", ColorFormat::HSV},
                    {"YUV444", ColorFormat::YUV444},
                    {"CFA", ColorFormat::CFA},
                    {"Unknown", ColorFormat::Unknown} };

                    return findValueInMap(formatMap, pixelLayout);
            }

            static std::string getColorFormatName(const images::datatypes::ColorFormat format) 
            {

                using namespace images::datatypes;

                typedef std::map<ColorFormat, std::string> ColorFormatNameMap;

                static const ColorFormatNameMap formatNameMap {
                    {ColorFormat::BGRA, "BGRA"},
                    {ColorFormat::RGB, "RGB"},
                    {ColorFormat::RGBA, "RGBA"},
                    {ColorFormat::Grayscale, "Grayscale"},
                    {ColorFormat::CMYK, "CMYK"},
                    {ColorFormat::Spectral, "Spectral"},
                    {ColorFormat::Lab, "Lab"},
                    {ColorFormat::HSV, "HSV"},
                    {ColorFormat::YUV444, "YUV444"},
                    {ColorFormat::CFA, "CFA"},
                    {ColorFormat::Unknown, "Unknown"} };

                    return findValueInMap(formatNameMap, format);
            }

            static uint8_t getColorFormatValue(const images::datatypes::ColorFormat format)
            {
                return static_cast<uint8_t>( format );
            }

            static images::datatypes::Layout getLayoutFromName(const std::string& layout) 
            {

                using namespace images::datatypes;
                
                typedef std::map<std::string, Layout> LayoutMap;

                static const LayoutMap layoutMap {
                    {"RowMajorInterleaved", Layout::RowMajorInterleaved},
                    {"ColumnMajorPlanar", Layout::ColumnMajorPlanar} };

                    return findValueInMap(layoutMap, layout);
            }

            static std::string getLayoutName(const images::datatypes::Layout layout) 
            {

                using namespace images::datatypes;

                typedef std::map<Layout, std::string> LayoutNameMap;

                static const LayoutNameMap layoutNameMap {
                    {Layout::RowMajorInterleaved, "RowMajorInterleaved"},
                    {Layout::ColumnMajorPlanar, "ColumnMajorPlanar"} };

                    return findValueInMap(layoutNameMap, layout);
            }

            static uint8_t getLayoutValue(const images::datatypes::Layout layout)
            {
                return static_cast<uint8_t>( layout );
            }                        
        }
    }
}


#endif  //imagesdatatypesmwimage_utils_header