#ifndef __FILE_CONTROL_HPP__
#define __FILE_CONTROL_HPP__

#include <string>
#include <direct.h>

class MakeFolder
{
public:

  MakeFolder(void){}

  virtual ~MakeFolder(void){}

  static bool checkExistenceOfFolder(const std::string folder_name) {
    if( _mkdir( folder_name.c_str() ) == 0 ){
      return true;
    } else {
      return false;
    }
  }
};


#endif // __FILE_CONTROL_H__