/**
* @file         ImageMetadata.hpp
*
* Purpose       Class definition of ImageMetadata class
*
* Authors       Dinesh Iyer
*
* Copyright     2019-2020 MathWorks, Inc.
*
*/

#ifndef imagesdatatypesimagemetadata_header
#define imagesdatatypesimagemetadata_header

#ifdef IMAGE_IN_CODEGEN_MODE
#include "datatypes_util.hpp" /* To deal with the fact that PackNGo has no include file hierarchy during test */
#include "ImageDefs.hpp"
#else
#include "images/datatypes/datatypes_util.hpp"
#include "images/datatypes/ImageDefs.hpp"
#endif

namespace images
{
    namespace datatypes
    {
        class ImageMetadata
        {
            public:
                // ------------------------------------------------
                // Life time
                // ------------------------------------------------
 
                /**
                  * Default constructor
                  * This represents metadata for an empty Image
                  */
                ImageMetadata()
                    : ImageMetadata( 0, 0, 3, 
                                     ClassUnderlying::Uint8,
                                     ColorFormat::RGB,
                                     Layout::RowMajorInterleaved )
                {}

                /**
                  * Create the metadata for describing an image
                  *
                  * @param - rows - Number of rows in the image
                  * @param - columns - Number of columns in the image
                  * @param - numChannels - Number of channels in the image
                  * @param - dtype - Datatype of the underlying data
                  * @param - format - ColorFormat of the image data
                  * @param - layout - Layout of the image data
                  * @param - relinquish - YES indicates the instance can give up its
                  *                        underlying data buffer without a data copy
                  */
                ImageMetadata( DimType rows,
                               DimType columns,
                               DimType numChannels,
                               ClassUnderlying dtype,
                               ColorFormat format,
                               Layout layout )
                    : fRows(rows), fColumns(columns), fNumChannels(numChannels),
                      fClassUnderlying(dtype),
                      fColorFormat(format),
                      fLayout(layout)
                {}

                // Not really necessary. But adding for completeness
                ~ImageMetadata() noexcept 
                {}

                // ------------------------------------------------
                // Copy and assignment operations
                // ------------------------------------------------
                // Miranda operations are sufficient for this class

                // ------------------------------------------------
                // Property Getters/Setters
                // ------------------------------------------------
                DimType getRows() const { return fRows; }
                DimType getColumns() const { return fColumns; }
                DimType getNumChannels() const { return fNumChannels; }
 
                ClassUnderlying getClassUnderlying() const { return fClassUnderlying; }
                ColorFormat getColorFormat() const { return fColorFormat; }
                Layout getLayout() const { return fLayout; }

                void setRows(DimType rows) { fRows = rows; }
                void setColumns(DimType cols) { fColumns = cols; }
                void setNumChannels(DimType chans) { fNumChannels = chans; }

                void setClassUnderlying(ClassUnderlying c) { fClassUnderlying = c; }
                void setColorFormat(ColorFormat cf) { fColorFormat = cf; }
                void setLayout(Layout l) { fLayout = l; }


            private:
                // ------------------------------------------------
                // Properties
                // ------------------------------------------------
                DimType fRows;
                DimType fColumns;
                DimType fNumChannels;
                ClassUnderlying fClassUnderlying;
                ColorFormat fColorFormat;
                Layout fLayout;
        };
    }
}


#endif // imagesdatatypesimagemetadata_header