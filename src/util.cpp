#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include "util.h"

std::string util::exec(char const* cmd) {
  char buffer[128];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 128, pipe) != NULL)
	result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return result;
}


const std::string util::revComp (const std::string & sequence){
  std::string newString = "";
  //cout << "Start - " << Sequence << "\n";
  for(int i = sequence.size()-1; i>=0; i+= -1) {
    char C = sequence.c_str()[i];
    if (C == 'A')
      {newString += 'T';}
    else if (C == 'C')
      {newString += 'G';}
    else if (C == 'G')
      {newString += 'C';}
    else if (C == 'T')
      {newString += 'A';}
    else if (C == 'N')
      {newString += 'N';}
    else {std::cout << "\nERROR IN RevComp - " << C << "\n";}
  }
  return newString;
}

const std::string util::baseName(const std::string & path){
  return path.substr(path.find_last_of("/\\")+1);
}
