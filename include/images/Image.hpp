/**
* @file         Image.hpp
*
* Purpose       Class declaration of Image class
*
* Authors       Dinesh Iyer
*
* Copyright     2019-2020 MathWorks, Inc.
*
*/

#ifndef imagesdatatypesimage_header
#define imagesdatatypesimage_header

#ifdef IMAGE_IN_CODEGEN_MODE
#include "datatypes_util.hpp" /* To deal with the fact that PackNGo has no include file hierarchy during test */
#include "ImageDefs.hpp"
#include "ImageMetadata.hpp"
#else
#include "images/datatypes/datatypes_util.hpp"
#include "images/datatypes/ImageDefs.hpp"
#include "images/datatypes/ImageMetadata.hpp"
#endif

#include <memory>

namespace images
{
    namespace datatypes
    {
        class IMAGESDATATYPESIMAGE_API Image
        {
            private:
                // ------------------------------------------------
                // Typedef's used
                // ------------------------------------------------
                typedef std::default_delete<BufferType[]> RawDataDeleter;
                
                // ------------------------------------------------
                // Interface of class that manages the lifetime of
                // the underlying data
                // ------------------------------------------------
                class IDataUnderlyingMgr
                {
                    public:
                        virtual ~IDataUnderlyingMgr() noexcept {}
                         
                        // Gain access to the raw data buffer
                        virtual BufferType* getDataUnderlying() const = 0;
 
                        // Give-up ownership of the underlying data buffer
                        virtual BufferType* releaseDataUnderlying() = 0;
                };
                
                // ------------------------------------------------
                // Takes ownership of the underling data.
                // Memory for underlying data is de-allocated if no
                // references to it remain.
                // ------------------------------------------------
                class DataUnderlyingMgrOwn : public IDataUnderlyingMgr
                {
                    public:
                        DataUnderlyingMgrOwn(BufferType* src) 
                            :fDataUnderlying(nullptr) 
                        {
                            if( src == nullptr ) 
                            {
                                return;
                            }

                            fDataUnderlying.reset(src, RawDataDeleter());
                        }

                        ~DataUnderlyingMgrOwn() noexcept {}

                        BufferType* getDataUnderlying() const override { return fDataUnderlying.get(); }
                        
                        BufferType* releaseDataUnderlying() override 
                        {
                            //  @TODO
                            return nullptr;
                        }

                    private:
                        std::shared_ptr<BufferType> fDataUnderlying;
                };
                
                // ------------------------------------------------
                // Uses the pointer provided as the underlying data
                // Does not dellocate memory.
                // ------------------------------------------------
                class DataUnderlyingMgrUseOnly : public IDataUnderlyingMgr
                {
                    public:
                        DataUnderlyingMgrUseOnly(BufferType* src)
                            :fDataUnderlying(src) 
                        {}

                        ~DataUnderlyingMgrUseOnly() noexcept {}

                        BufferType* getDataUnderlying() const override { return fDataUnderlying; }

                        BufferType* releaseDataUnderlying() override 
                        { 
                            auto tempBuffer = fDataUnderlying; 
                            fDataUnderlying = nullptr;
                            return tempBuffer;
                        }

                    private:
                        BufferType* fDataUnderlying;
                };
                
                // ------------------------------------------------
                // Takes ownership of the underling data.
                // Supports relinquishing the underlying data without
                // doing a deep copy of the data
                // ------------------------------------------------
                class DataUnderlyingMgrSingleUse : public IDataUnderlyingMgr
                {
                    public:
                        DataUnderlyingMgrSingleUse(BufferType* src)
                            : fDataUnderlying( src, RawDataDeleter() )
                        {}

                        ~DataUnderlyingMgrSingleUse() noexcept {}

                        BufferType* getDataUnderlying() const override { return fDataUnderlying.get(); }

                        BufferType* releaseDataUnderlying() override { return fDataUnderlying.release(); }

                    private:
                        std::unique_ptr<BufferType, RawDataDeleter> fDataUnderlying;
                };
            
            public:
                // ------------------------------------------------
                // Life time
                // ------------------------------------------------
 
                /**
                  * Default constructor
                  * This represent an empty image
                  * Also required by MATLAB Coder
                  */
                Image();     
                
                /**
                  * Construct an image.
                  *
                  * Allow users to create an image with certain optimizations that
                  * can avoid unnecessary copies of the underlying data.
                  *
                  * The image created with RelinquishDataUnderlying::YES can be copied and assigned.
                  * However, when copied or assigned, the target image will contain a deep copy
                  * of the data of the source image. For IM1 created with this flag = TRUE
                  * Image im2 = im1; % IM2 will have its own copy of the data buffer i.e
                  * deep copy will be done
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
                Image( DimType rows,
                       DimType columns,
                       DimType numChannels,
                       ClassUnderlying dtype,
                       ColorFormat format,
                       Layout layout,
                       RelinquishDataUnderlying relinquish = RelinquishDataUnderlying::NO );

                /**
                  * Construct an image.
                  *
                  * Allow users to create an image with certain optimizations that
                  * can avoid unnecessary copies of the underlying data.
                  *
                  * The image created with RelinquishDataUnderlying::YES can be copied and assigned.
                  * However, when copied or assigned, the target image will contain a deep copy
                  * of the data of the source image. For IM1 created with this flag = TRUE
                  * Image im2 = im1; % IM2 will have its own copy of the data buffer i.e
                  * deep copy will be done
                  *
                  * @param - metadata - Instance of ImageMetadata class 
                  * @param - relinquish - YES indicates the instance can give up its
                  *                        underlying data buffer without a data copy
                  */
                Image( const ImageMetadata& metadata,
                       RelinquishDataUnderlying relinquish = RelinquishDataUnderlying::NO );
                
                /**
                  * Construct an image.
                  *
                  * Allow users to create an image by specifying the dimensions, datatype,
                  * memory arrangement and also provide the buffer containing the image
                  * data
                  *
                  *
                  * @param - src - Source buffer containing data. Assumed that src has
                  *                sufficient number of bytes to store the entire image
                  *                data
                  * @param - rows - Number of rows in the image
                  * @param - columns - Number of columns in the image
                  * @param - numChannels - Number of channels in the image
                  * @param - dtype - Datatype of the underlying data
                  * @param - format - ColorFormat of the image data
                  * @param - layout - Layout of the image data
                  * @param - mgmt - Determine how the input source of data needs to be
                  *                 handled by the image. By default, a deep copy of the
                  *                 source is made.
                  * @param - relinquish - YES indicates the instance can give up its
                  *                        underlying data buffer without a data copy
                  */
                Image( void* src, 
                       DimType rows,
                       DimType columns,
                       DimType numChannels,
                       ClassUnderlying dtype,
                       ColorFormat format,
                       Layout layout,
                       SourceMgmt mgmt = SourceMgmt::COPY_SRC,
                       RelinquishDataUnderlying relinquish = RelinquishDataUnderlying::NO );

                /**
                  * Construct an image.
                  *
                  * Allow users to create an image by specifying the dimensions, datatype,
                  * memory arrangement and also provide the buffer containing the image
                  * data
                  *
                  *
                  * @param - src - Source buffer containing data. Assumed that src has
                  *                sufficient number of bytes to store the entire image
                  *                data
                  * @param - rows - Number of rows in the image
                  * @param - columns - Number of columns in the image
                  * @param - numChannels - Number of channels in the image
                  * @param - dtype - Datatype of the underlying data
                  * @param - format - ColorFormat of the image data
                  * @param - layout - Layout of the image data
                  */
                Image( void const* dataUnderlying,
                       DimType rows,
                       DimType columns,
                       DimType numChannels,
                       ClassUnderlying dtype,
                       ColorFormat format,
                       Layout layout );

                /**
                  * Construct an image.
                  *
                  * Allow users to create an image by specifying the dimensions, datatype,
                  * memory arrangement and also provide the buffer containing the image
                  * data
                  *
                  *
                  * @param - src - Source buffer containing data. Assumed that src has
                  *                sufficient number of bytes to store the entire image
                  *                data
                  * @param - metadata - Instance of ImageMetadata class 
                  * @param - mgmt - Determine how the input source of data needs to be
                  *                 handled by the image. By default, a deep copy of the
                  *                 source is made.
                  * @param - relinquish - YES indicates the instance can give up its
                  *                        underlying data buffer without a data copy
                  */
                Image( void* src, 
                       const ImageMetadata& metadata,
                       SourceMgmt mgmt = SourceMgmt::COPY_SRC,
                       RelinquishDataUnderlying relinquish = RelinquishDataUnderlying::NO );
                       
                // Not really necessary. But adding for completeness
                ~Image() noexcept;
                
                // ------------------------------------------------
                // Copy and assignment operations
                // ------------------------------------------------
 
                /**
                  * Perform a shallow copy of the image
                  */
                Image(const Image& otherImage);

                Image& operator=(const Image& otherImage);
 
                // Move operations
                Image(Image&& otherImage) noexcept;
                
                Image& operator=(Image&& otherImage) noexcept;
                
                // ------------------------------------------------
                // Public methods
                // ------------------------------------------------
                /**
                  * Creates a deep copy i.e. allocates a new memory buffer
                  * and copies data into it.
                  *
                  * @param - relinquish - YES indicates the instance can give up its
                  *                        underlying data buffer without a data copy
                  * @return - Image
                  */
                Image deepCopy(RelinquishDataUnderlying relinquish = 
                                            RelinquishDataUnderlying::NO) const;
 
                /**
                  * Access the underlying data buffer
                  *
                  * Using the data buffer, the user can modify the contents of the data.
                  * As the internal buffer is a shared resource, changing it can have
                  * indeterminate consequences.
                  * So, if the user wants to modify the buffer, the recommendation is
                  * to create a deepCopy() FIRST.
                  *
                  * @param - none
                  * @return - Pointer of type BufferType
                  */
                BufferType* getDataUnderlying() const;

                BufferType* getUniqueDataUnderlying();

                /**
                  * Query whether calling releaseDataUnderlying() makes a deep copy
                  * of the underlying data before passing it to the user
                  *
                  * This is necessary if the user wants to track when a memcpy is being
                  * done when working with Image instances. The table specified in this
                  * spec contains details.
                  *
                  * @param - none
                  * @return - bool. TRUE if a deep copy will be done when calling releaseDataUnderlying
                  */
                bool isDeepCopyOfDataOnRelease() const;
     
                /**
                  * Give up ownership of the underlying data buffer
                  *
                  * This is useful if the user wants to pass the data to another
                  * object that can take ownership of data buffers.
 
                  * NOTE: Using the Image instance to access image data after calling
                  *       this method have indeterminate behaviour.
                  *
                  * @param - none
                  * @return - Pointer of type BufferType
                  */
                BufferType* releaseDataUnderlying();

                /**
                  * Resize the image
                  *
                  * This method "blows out" any valid image data that is held by the 
                  * image instance and allocates a new buffer that has sufficient
                  * space to hold image data for the newer dimensions.
                  * This method allows changing only the number of rows and/or columns
                  * AND not the number of channels or any other properties.
                  *
                  * @param - newNumRows - Number of rows in the target image
                  * @param - newNumCols - Number of columns in the target image
                  * @return - none
                  */
                void resize(DimType newNumRows, DimType newNumCols);

                size_t computeNumBytesInImage() const;
                
                // ------------------------------------------------
                // Property Getters/Setters
                // ------------------------------------------------
                DimType getRows() const { return fMetadata.getRows(); }
                DimType getColumns() const { return fMetadata.getColumns(); }
                DimType getNumChannels() const { return fMetadata.getNumChannels(); }
 
                ClassUnderlying getClassUnderlying() const { return fMetadata.getClassUnderlying(); }
                ColorFormat getColorFormat() const { return fMetadata.getColorFormat(); }
                Layout getLayout() const { return fMetadata.getLayout(); }
                
            private:
                // ------------------------------------------------
                // Private Helper functions
                // ------------------------------------------------
                void copyProps(const Image& otherImage);
                
                BufferType* allocateBufferForDataUnderlying() const;
                
                IDataUnderlyingMgr* createDataUnderlyingMgr( BufferType* src, 
                                                             SourceMgmt mgmt = SourceMgmt::COPY_SRC,
                                                             RelinquishDataUnderlying relinquish = 
                                                                            RelinquishDataUnderlying::NO );

                std::unique_ptr<BufferType, RawDataDeleter> copyRawData(const BufferType* src) const;
                
                // ------------------------------------------------
                // Properties
                // ------------------------------------------------                
                ImageMetadata fMetadata;

                // Instance of class that manages the raw image data
                std::shared_ptr<IDataUnderlyingMgr> fDataUnderlyingMgr;
 
                // Tracks whether the instance can relinquish underlying data
                // without performing a deep copy.
                RelinquishDataUnderlying fRelinquish;
        };

        inline Image makeImage(void const* dataUnderlying,
                               DimType rows,
                               DimType columns,
                               DimType numChannels,
                               ClassUnderlying classUnderlying,
                               ColorFormat colorFormat,
                               Layout layout) {
            return Image(
                dataUnderlying,
                rows,
                columns,
                numChannels,
                classUnderlying,
                colorFormat,
                layout);
        }

        inline Image makeImage(DimType rows,
                               DimType columns,
                               DimType numChannels,
                               ClassUnderlying classUnderlying,
                               ColorFormat colorFormat,
                               Layout layout) {
            return Image(
                rows,
                columns,
                numChannels,
                classUnderlying,
                colorFormat,
                layout);
        }
    }
}

#endif // imagesDatatypesimage_header