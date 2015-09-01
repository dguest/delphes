#ifndef H5_TYPES_HH
#define H5_TYPES_HH

#include "H5Cpp.h"
#include <vector>
#include <string>

#define H5_INSERT(INTO, INTO_CLASS, MEMBER)				\
  do {									\
    using h5::insert;							\
    INTO_CLASS dummy;							\
    insert(INTO, #MEMBER, offsetof(INTO_CLASS, MEMBER), dummy.MEMBER);	\
  } while (false)

namespace h5 {
  H5::DataType type(int);
  H5::DataType type(double);
  template <typename T>
  H5::DataType type(const h5::vector<T>& dummy) {
    T subdummy;
    const auto subtype = type(subdummy);
    return H5::VarLenType(&subtype);
  }

  template <typename T>
  void insert(H5::CompType& into, const std::string& name,
	      size_t offset, T dummy) {
    into.insertMember(name, offset, type(dummy));
  }

}

#endif
