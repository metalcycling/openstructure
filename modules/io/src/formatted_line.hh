#ifndef OST_IO_FORMATTED_LINE_HH
#define OST_IO_FORMATTED_LINE_HH

#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#include <cstdio>

#include <ost/string_ref.hh>

/*
  Author: Marco Biasini
 */

namespace ost { namespace io {

namespace fmt {

struct LPadded : public StringRef {
  LPadded(const char* s, size_t l): 
    StringRef(s, l)
  { }
  
  explicit LPadded(const String& s): StringRef(s.data(), s.size())
  { }
  
  explicit LPadded(const char* s): StringRef(s, strlen(s))
  { }
};

struct RPadded : public StringRef {
  RPadded(const char* s, size_t l): 
    StringRef(s, l)
  { }
  
  explicit RPadded(const String& s): StringRef(s.data(), s.size())
  { }
  
  explicit RPadded(const char* s): StringRef(s, strlen(s))
  { }
};


struct LPaddedInt {
  LPaddedInt(int val) {
    size_t curr=0;
    bool minus=val<0;
    val=std::abs(val);
    if (minus) {
      data[curr]='-';
      ++curr;
    }
    do {
      int m=val%10;
      data[curr]='0'+m;
      ++curr;
      val/=10;
    } while(val);
    // swap
    for (size_t i=0, e=(curr-int(minus))/2; i<e; ++i) {
      std::swap(data[int(minus)+i], data[curr-i-1]);
    }
    data[curr]='\0';
    len=curr;
  }
  char   data[20];
  size_t len;
};

// BIG FAT WARNING: Before using this class to output floats with lots of 
// digits, make sure you indeed get the result you want. The function has 
// not been tested enough for numerical stability.
struct LPaddedFloat {
  LPaddedFloat(Real val, int prec)
  {
    switch(prec){
      case 0:
        len = snprintf(data, sizeof(data), "%.0f", val);
        break;
      case 1:
        len = snprintf(data, sizeof(data), "%.1f", val);
        break;
      case 2:
        len = snprintf(data, sizeof(data), "%.2f", val);
        break;
      case 3:
        len = snprintf(data, sizeof(data), "%.3f", val);
        break;
      default:
        Real rounded_val=round(val*pow(Real(10), prec))*pow(Real(0.1), prec);
        size_t curr=0;
        bool minus=rounded_val<0;
        rounded_val=std::abs(rounded_val);
        int int_val=int(rounded_val);
        if (minus) {
          data[curr]='-';
          ++curr;
        }
        do {
          int m=int_val%10;
          data[curr]='0'+m;
          ++curr;
          int_val/=10;
        } while(int_val);
        // swap
        for (size_t i=0, e=(curr-int(minus))/2; i<e; ++i) {
          std::swap(data[int(minus)+i], data[curr-i-1]);
        }
        data[curr]='\0';
        len=curr;
        if (prec==0) {
          return;
        }
        data[curr]='.';
        curr++;
        rounded_val-=int(rounded_val);
        while(prec>0) {
          rounded_val*=10;
          int m=int(rounded_val);
          rounded_val-=int(rounded_val);
          data[curr]='0'+m;
          curr++;
          --prec;
        }
        data[curr]='\0';
        len=curr;
    }
  }
  char   data[20];
  size_t len;  
};



}

class LineSlice {
public:
  LineSlice(char* data, size_t l): data_(data), len_(l) 
  {
  }

  LineSlice& operator=(const StringRef& str)
  {
    assert(str.length()==len_);
    memcpy(data_, str.data(), str.size());
    return *this;
  }
  
  LineSlice& operator=(const fmt::LPadded& str)
  {
    assert(str.size()<=len_);
    memcpy(data_+len_-str.size(), str.data(), str.size());
    return *this;
  }
  
  LineSlice& operator=(const fmt::RPadded& str)
  {
    assert(str.size()<=len_);
    memcpy(data_, str.data(), str.size());
    return *this;
  }
  
  LineSlice& operator=(const fmt::LPaddedInt& i) 
  {
    assert(i.len<=len_);
    memcpy(data_+len_-i.len, i.data, i.len);
    return *this;
  }
  
  LineSlice& operator=(const fmt::LPaddedFloat& f)
  {
    assert(f.len<=len_);
    memcpy(data_+len_-f.len, f.data, f.len);
    return *this;
  }
  void Clear()
  {
    memset(data_, ' ', len_);
  }
private:
  char* data_;
  size_t len_;
};

class FormattedLine {
public:
  FormattedLine(size_t width): 
    data_(new char[width]), len_(width)
  {
    this->Clear();
  }
  
  void Clear()
  {
    memset(data_, ' ', len_);
  }
  ~FormattedLine() { delete[] data_; }
  
  LineSlice operator()(int start, int len) 
  { 
    assert(start>=0 && start+len<=static_cast<int>(len_));
    return LineSlice(data_+start, len);
  }
  
  const char* Data() const { return data_; }
  
  size_t GetWidth() const { return len_; }
  
  char operator[](size_t index) const
  {
    assert(index<len_);
    return data_[index];
  }
  
  char& operator[](size_t index)
  {
    assert(index<len_);
    return data_[index];
  }
private:
  FormattedLine& operator=(const FormattedLine& rhs);
  FormattedLine(const FormattedLine& rhs);
  char*       data_;
  size_t      len_;
};

inline std::ostream& operator<<(std::ostream& stream, const FormattedLine& line)
{
  stream.write(line.Data(), line.GetWidth()); 
  stream << std::endl;
  return stream;
}


}}

#endif
