#include <iostream>

#include <ost/log.hh>
#include <ost/img/alg/stat.hh>
#include "highest_peak_search_3d.hh"

namespace ost { namespace img { namespace alg {

namespace detail {

  bool is_within( const Point& origin, const Point& p, int dist)
  {
    int dist2=dist*dist;
    int dx=origin[0]-p[0];
    if(dx*dx<dist2) {
      int dy=origin[1]-p[1];
      if(dy*dy<dist2) {
        int dz=origin[2]-p[2];
        if((dx*dx+dy*dy+dz*dz)<dist2) {
          return true;
        }
      }
    }
    return false;
  }

  struct is_smaller_and_within_radius
  {
    is_smaller_and_within_radius(const Peak& input_cand_peak, int input_exclusion_radius): cand_peak(input_cand_peak),
      exclusion_radius(input_exclusion_radius) {}

    bool operator()(Peak& peak) const
    {
      return (is_within(cand_peak,peak,exclusion_radius) && cand_peak.GetValue() > peak.GetValue());
    }

    Peak cand_peak;
    int exclusion_radius;
  };

  //! Collects and manages a list of the highest peaks in a 3D correlation map
  class PeakCollector
  {

    public:

      //! Constructs the collector and initializes it with user defined parameters
      PeakCollector(int max_num_peaks, int exclusion_radius):
        max_num_peaks_(max_num_peaks),
        exclusion_radius_(exclusion_radius)
      {}

      //! Updates peak list adding the candidate peak if the conditions are met
      /*!
        given a new peak, it will be added if it is outside
        exclusion_radius of any other existing peak
        if it is inside exclusion_radius and its value is larger
        that any other ones inside, it will be added, and all smaller
        peaks within exclusion_radius will be removed
        from the internal list. Only if there is a larger peak within
        exclusion_radius will this one be ignored.

        the decision wether a peak is at all to be considered is in
        the hands of the caller

        if a peak gets added, and the list has reached max_size, the
        peak with the lowest value will be ejected

        returns true is the peak has been added, false otherwise
      */
      bool Add(const img::Peak& cand_peak)
      {
        static img::Peak minpeak(img::Point(0,0,0),-std::numeric_limits<Real>::max());
        bool insert=true;
        bool remove_smaller=false;
        img::PeakList::iterator min_it;

        if(!list_.empty()) {
          min_it = list_.begin();

          for(img::PeakList::iterator it=list_.begin();it!=list_.end();++it)
          {
            if(it->GetValue() < min_it->GetValue()) min_it = it;
            if (it->GetValue() >= cand_peak.GetValue() && detail::is_within(cand_peak, *it, exclusion_radius_)) {
              insert=false;
            } else if (cand_peak.GetValue() >= it->GetValue() && detail::is_within(cand_peak, *it, exclusion_radius_)) {
              remove_smaller=true;
            }
          }

          if (!insert) return false;

          if (cand_peak.GetValue() < min_it->GetValue()) {
            return false;
          }

          if (remove_smaller)
          {
            img::PeakList::iterator new_end = remove_if(list_.begin(),list_.end(),
            detail::is_smaller_and_within_radius(cand_peak,exclusion_radius_));
            list_.erase(new_end,list_.end());
          }

          if (list_.size() == max_num_peaks_) {
            list_.erase(min_it);
          }
        }

        list_.push_back(cand_peak);
        return true;
      }

       //! Gets list of peaks from collector
      img::PeakList GetPeakList() const { return list_; }

      //! Clears the peak list removing all peaks in the list
      void ClearPeakList() { list_.clear(); }

    private:

      img::PeakList list_;      
      unsigned int max_num_peaks_;
      int exclusion_radius_;

  };

} // ns

template <typename T, class D>
  void HighestPeakSearch3DBase::VisitState(const ImageStateImpl<T,D>& isi)
  {

    LOG_DEBUG("highest peak search (" 
              << "max number of peaks: " << max_num_peaks_ << ", "
              << "exclusion radius: " << exclusion_radius_ << ", "
              << "threshold: " << threshold_ << ")");

    if (max_num_peaks_ < 1)
    {
      throw Error("Error: Maximum number of peaks smaller than 1. Cannot create peak list");
    }

    if(ext_list_.size()==0)
    {
      LOG_DEBUG(" excluded regions: none");
    } 
    else 
    {
      std::stringstream ss;
      for(ExtList::const_iterator it=ext_list_.begin();it!=ext_list_.end();++it) 
      {
        ss << " " << (*it);
      }
      LOG_DEBUG("excluded regions:" << ss.str());
    }

    detail::PeakCollector peak_collector(max_num_peaks_,exclusion_radius_);

    Extent img_ext = isi.GetExtent();
    for(ExtentIterator img_it(img_ext);!img_it.AtEnd();++img_it) 
    {
      if(is_excluded(Point(img_it))) continue;
    
      Real val = Val2Val<T,Real>(isi.Value(img_it));
      if( val >= threshold_ )
      {
        Peak peak(Point(img_it),val);
        peak_collector.Add(peak);
      } 
  
    }
    peaks_=peak_collector.GetPeakList();
    peak_collector.ClearPeakList();
  }


} // ns


template class TEMPLATE_DEF_EXPORT image_state::ImageStateNonModAlgorithm<alg::HighestPeakSearch3DBase>;


}}  // ns
