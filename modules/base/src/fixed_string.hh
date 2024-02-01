#ifndef OST_SMALL_STRING_HH
#define OST_SMALL_STRING_HH

#include <ost/base.hh>
#include <cassert>
#include <string.h>

namespace ost {

/// \brief string class that uses an array of static size to hold the 
///        characters
template <size_t LENGTH>
class TEMPLATE_EXPORT FixedString {
public:
  FixedString() {
    ::memset(bytes_, 0, LENGTH+1);
  }

  explicit FixedString(const String& str) {
    this->assign(str);
  }

  explicit FixedString(const char* str) {
    this->assign(str);
  }


  size_t length() const {
    return strlen(bytes_);
  }

  size_t size() const {
    return this->length();
  }

  char operator[](size_t index) const {
    return bytes_[index];
  }

  char& operator[](size_t index) {
    return bytes_[index];
  }

  bool operator==(const String& rhs) const {
    return  !strcmp(bytes_, rhs.c_str());
  }

  bool operator!=(const String& rhs) const {
    return !this->operator==(rhs);
  }
  bool operator==(const FixedString<LENGTH>& rhs) const {
    return  !strcmp(bytes_, rhs.bytes_);
  }

  bool operator!=(const FixedString<LENGTH>& rhs) const {
    return !this->operator==(rhs);
  }
  size_t capacity() const {
    return LENGTH;
  }
  FixedString<LENGTH>& operator=(const String& rhs) {
    this->assign(rhs);
    return *this;
  }

  template <typename DS>
  void Serialize(DS& ds) {
    RawSerialize(ds, bytes_, LENGTH);
  }
  const char* c_str() const {
    return bytes_;
  }
  char* data() { return bytes_; }
  const char* data() const { return bytes_; }
private:
  void assign(const String& str) {
    assert(str.length()<=LENGTH);
    strcpy(bytes_, str.c_str());
  }

  void assign(const char* str) {
    assert(strlen(str)<=LENGTH);
    strcpy(bytes_, str);
  }
  char bytes_[LENGTH+1];
};

}

#endif
