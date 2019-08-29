#ifndef __itkFastGrowCutSegmentationImageFilter_h
#define __itkFastGrowCutSegmentationImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkSimpleDataObjectDecorator.h"
#include "itkVectorContainer.h"
//#include "itkCommand.h"

#include "vcl_list.h"

#include <vcl_vector.h>

#include "FastGrowCutSegmenter.h"

//#include "itkGrowCutSegmentationUpdateFilter.h"

namespace itk
{

/**  \class GrowCutSegmentationFilter
 * \brief Segmentation or multiple objects based on a set of gestures.
 *
 * Given a mask image containing the gestures for foreground classes
 * and their background, employs grow cut segmentation to produce
 * segmentation. Supports passing an existing segmentation with
 * additional gestures for editing a segmentation produced as a
 * result of the same or different algorithm.
 *
 * The filter is based on the paper "GrowCut:Interactive Multi-Label
 * N-D Image Segmentation By Cellular Automata", Vladimir Vezhnevets,
 * Vadim Konouchine
 *
 * Modified Version: The inputs consist of the input intensity image,
 * optionally the input Label Image, and the Seed image.  The Label
 * image and the Seed image default to using unsigned chars and
 * double respectively.
 *
 * This algorithm is implemented scalar images. Vector Images are not
 * supported.
 *
 *
**/


/* template<class TInputImage,  */
/*   class TOutputImage, class TLabelPixelType = short,  */
/*   class TSeedPixelType = float >  */
template<class TInputImage,
         class TOutputImage>
class FastGrowCutSegmentationImageFilter: public ImageToImageFilter<TInputImage, TOutputImage>
{

public:
    /** Standard class type alias. */
    using Self = FastGrowCutSegmentationImageFilter;
    using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
    using Pointer = SmartPointer<Self>;
    using ConstPointer = SmartPointer<const Self>;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods).  */
    itkTypeMacro(FastGrowCutSegmentationImageFilter,
                 ImageToImageFilter);

    /** Image related type alias. */
    itkStaticConstMacro(ImageDimension, unsigned int,
                        TInputImage::ImageDimension );

    using InputImageType = TInputImage;
    using InputImagePointer = typename InputImageType::Pointer;
    using InputImageConstPointer = typename InputImageType::ConstPointer;

    using InputPixelType = typename InputImageType::PixelType;
    using InputIndexType = typename InputImageType::IndexType;
    using SizeType = typename InputImageType::SizeType;

    using OutputImageType = TOutputImage;
    using OutputImagePointer = typename OutputImageType::Pointer;
    using OutputImageRegionType = typename OutputImageType::RegionType;
    using OutputPixelType = typename OutputImageType::PixelType;
    using OutputIndexType = typename OutputImageType::IndexType;
    using OutputSizeType = typename InputImageType::SizeType;


    /** Smart Pointer type to a DataObject. */
    using DataObjectPointer = typename DataObject::Pointer;


    /** Image dimension constants */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);

    /** Index type alias support. */
    using IndexType = Index<itkGetStaticConstMacro(InputImageDimension)>;

    /** InputSizeType type alias support **/
    using InputSizeType = typename InputImageType::SizeType;

    /** SeedImage type alias support */
    /* indicates the Seed of a label for a given cell */
    using SeedImageType = TOutputImage;

    /** SeedImagePointer type alias support. */
    using SeedImagePointer = typename SeedImageType::Pointer;

    /** Set the Input Image **/
    void SetInputImage( const InputImageType *in)
    {
        this->ProcessObject::SetNthInput(0, const_cast< InputImageType *>(in) );
    }

    const InputImagePointer GetInputImage( );

    /** Set/Get the Label Image **/
    void SetLabelImage( const OutputImageType *f)
    {
        this->ProcessObject::SetNthInput(2, const_cast< OutputImageType *>(f) );
        m_LabelImage = static_cast< OutputImageType *>(this->ProcessObject::GetInput(2));
    }

    const OutputImagePointer GetLabelImage();

    /** Get the Seed image **/
    const SeedImagePointer GetSeedImage();

    /** Set the Seed Image **/
    void SetSeedImage( const SeedImageType *s )
    {
        this->ProcessObject::SetNthInput(1, const_cast< SeedImageType *>(s));
        //m_SeedImage = static_cast< SeedImageType *>(this->ProcessObject::GetInput(2));
    }

    const SeedImagePointer GetUpdatedSeedImage();

    itkSetMacro(InitializationFlag, bool);

    //processing functions
    void Initialization();
    void RunFGC();

protected:

    FastGrowCutSegmentationImageFilter();
    ~FastGrowCutSegmentationImageFilter() {};

    // Override since the filter needs all the data for the algorithm
    void GenerateInputRequestedRegion();

    // Override since the filter produces the entire dataset
    void EnlargeOutputRequestedRegion(DataObject *output);

    void GenerateData();

    void Initialize(OutputImageType *output);

    void PrintSelf ( std::ostream &os, Indent indent ) const;

    void GrowCutSlowROI( TOutputImage *);


private:

    FastGrowCutSegmentationImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &); // purposely not implemented

    void GetRegionOfInterest();

    void ComputeLabelVolumes(TOutputImage *outputImage, vcl_vector< unsigned > &volumes, vcl_vector< unsigned > &phyVolumes);

    OutputImagePointer                       m_LabelImage;
    SeedImagePointer                         m_SeedImage;


    std::vector<OutputPixelType> m_imSeedVec;
    std::vector<OutputPixelType> m_imLabVec;
    std::vector<InputPixelType> m_imSrcVec;
    std::vector<long> m_imROI;

    //logic code
    FGC::FastGrowCut<InputPixelType, OutputPixelType> *m_fastGC;

    //state variables
    bool m_InitializationFlag;

};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastGrowCutSegmentationImageFilter.hxx"
#endif

#endif
