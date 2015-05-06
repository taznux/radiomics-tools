#ifndef FAST_MARCHING_H
#define FAST_MARCHING_H

#include <math.h>
#include <queue>
#include <set>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iterator>

#include "utilities.h"

#include "fibheap.h"

namespace FGC
{

/***************************************************************************
* class HeapNode
***************************************************************************/
class HeapNode : public FibHeapNode
{
    float   N;
    long IndexV;

public:
    HeapNode() : FibHeapNode()
    {
        N = 0;
    }
    virtual void operator =(FibHeapNode &RHS);
    virtual int  operator ==(FibHeapNode &RHS);
    virtual int  operator <(FibHeapNode &RHS);
    virtual void operator =(double NewKeyVal );
    virtual void Print();
    double GetKeyValue()
    {
        return N;
    }
    void SetKeyValue(double n)
    {
        N = n;
    }
    long int GetIndexValue()
    {
        return IndexV;
    }
    void SetIndexValue( long int v)
    {
        IndexV = v;
    }
};

const float  DIST_INF = std::numeric_limits<float>::max();
const float  DIST_EPSION = 1e-3;
unsigned char NNGBH = 26;
typedef float FPixelType;

template<typename SrcPixelType, typename LabPixelType>
class FastGrowCut
{
public:
    FastGrowCut();
    ~FastGrowCut();

    void SetSourceImage(const std::vector<SrcPixelType> &imSrc);
    void SetSeedlImage(std::vector<LabPixelType> &imSeed);
    void SetWorkMode(bool bSegInitialized = false);
    void SetImageSize(const std::vector<long> &imSize);
    void DoSegmentation();
    void GetLabeImage(std::vector<LabPixelType> &imLab);
    void GetForegroundmage(std::vector<LabPixelType> &imFgrd);

private:
    void InitializationAHP();
    void DijkstraBasedClassificationAHP();


    std::vector<SrcPixelType> m_imSrc;
    std::vector<LabPixelType> m_imSeed;
    std::vector<LabPixelType> m_imLabPre;
    std::vector<FPixelType> m_imDistPre;
    std::vector<LabPixelType> m_imLab;
    std::vector<FPixelType> m_imDist;

    std::vector<long> m_imSize;
    long m_DIMX, m_DIMY, m_DIMZ, m_DIMXY, m_DIMXYZ;
    std::vector<int> m_indOff;
    std::vector<unsigned char>  m_NBSIZE;

    FibHeap *m_heap;
    HeapNode *m_hpNodes;
    bool m_bSegInitialized;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "FastGrowCutSegmenter.hxx"
#endif

#endif
