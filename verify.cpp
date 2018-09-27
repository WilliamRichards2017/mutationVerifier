#include <stdint.h>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

struct region{
  std::string startChrom;
  int32_t startPos;
  std::string endChrom;
  int32_t endPos;
};

int main(int argc, char *argv[]){
  std::cout << "inside main function bb " << std::endl;

  std::string regionString = argv[1];
 
  int16_t pos = regionString.find_first_of('-');
  std::string startR = regionString.substr(0, pos);
  std::string endR = regionString.substr(pos+1);

  int16_t startP = startR.find_first_of(':');
  int16_t endP = endR.find_first_of(':');

  
  region reg{.startChrom = startR.substr(0,startP), .startPos = std::stoi(startR.substr(startP+1)), .endChrom = endR.substr(0,endP), .endPos = std::stoi(endR.substr(endP+1))};
 

  std::cout << "region is: " << reg.startChrom << ':' << reg.startPos << '-' << reg.endChrom << ':' << reg.endPos << std::endl;

  return 0;

}
