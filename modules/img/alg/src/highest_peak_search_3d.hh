#ifndef OST_IMG_ALG_HIGHEST_PEAK_SEARCH_3D_H
#define OST_IMG_ALG_HIGHEST_PEAK_SEARCH_3D_H

#include <vector>

#include <ost/img/algorithm.hh>
#include <ost/img/peak.hh>
#include <ost/img/image_state.hh>
#include "module_config.hh"

namespace ost { namespace img { namespace alg {

class DLLEXPORT_IMG_ALG HighestPeakSearch3DBase
{
  public:

    //! Initialization
    /*!
     - number of peaks: number of peaks in the returned list
     - exclusion radius: radius of the exclusion volume around each peak (in voxels)
       (peaks will not be closer than this radius)
     - threshold: minimum value that a voxel must have to be considered a peak candidate
    */
        

    HighestPeakSearch3DBase(int max_num_peaks=10,
                            int exclusion_radius=10,
                            Real threshold=0):
      max_num_peaks_(max_num_peaks),
      exclusion_radius_(exclusion_radius),
      threshold_(threshold)
    {}
	         
    //! image state algorithm interface
    template <typename T, class D>
    void VisitState(const ImageStateImpl<T,D>& isi);
    
    static String GetAlgorithmName() {return "HighestPeakSearch3D";}
  
    //! Main "return" function returns a list of Peaks
    PeakList GetPeakList() const { return peaks_;}
  
    //! Clears the peak list removing all peaks in the list
    void ClearPeakList() { peaks_.clear(); }

  private:
  
    typedef std::vector<Extent> ExtList;

    int max_num_peaks_;
    int exclusion_radius_;
    Real threshold_;
    PeakList peaks_;
    ExtList ext_list_;

    bool is_excluded(const Point& p) 
    {
      for(ExtList::const_iterator it=ext_list_.begin();it!=ext_list_.end();++it) 
      {
        if(it->Contains(p)) return true;
      }
      return false;
    }
};

// injection of base into algorithmic class
typedef ImageStateNonModAlgorithm<HighestPeakSearch3DBase> HighestPeakSearch3D;

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateNonModAlgorithm<alg::HighestPeakSearch3DBase>)
}} // ns 

#endif // IPLT_ALG_DOCKING_HIGHEST_PEAK_SEARCH_3D_H
