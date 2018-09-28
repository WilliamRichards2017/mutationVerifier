#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <utility>

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

std::string util::baseName(std::string path){
  return path.substr(path.find_last_of("/\\")+1);
}
