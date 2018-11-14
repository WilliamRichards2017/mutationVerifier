#ifndef __SRC_UTIL__
#define __SRC_UTIL__

#include <string>

class util{
 public:
  static std::string exec(char const*);
  static const std::string baseName(const std::string &);
  static const std::string revComp(const std::string &);
  
};

#endif //__SRC_UTIL__
