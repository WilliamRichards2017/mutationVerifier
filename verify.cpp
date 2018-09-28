#include <stdint.h>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "util.h"

struct region{
  std::string startChrom;
  int32_t startPos;
  std::string endChrom;
  int32_t endPos;
};


void bamToFasta(std::string bamFile){
  std::string command = "samtools bam2fq ";
  command += bamFile;
  command += " | seqtk seq -A - > ";
  command += util::baseName(bamFile);
  command += ".fa";

  std::cout << "running command " << command << std::endl;

  util::exec(command.c_str());

}


int main(int argc, char *argv[]){

    std::string regionString = argv[1];
    std::string subjectBam = argv[2];
    
    int16_t pos = regionString.find_first_of('-');
    std::string startR = regionString.substr(0, pos);
    std::string endR = regionString.substr(pos+1);
    
    int16_t startP = startR.find_first_of(':');
    int16_t endP = endR.find_first_of(':');
    
    
    region reg{.startChrom = startR.substr(0,startP), .startPos = std::stoi(startR.substr(startP+1)), .endChrom = endR.substr(0,endP), .endPos = std::stoi(endR.substr(endP+1))};
    
    std::cout << "region is: " << reg.startChrom << ':' << reg.startPos << '-' << reg.endChrom << ':' << reg.endPos << std::endl;
    std::cout << "subject bam file is: " << subjectBam << std::endl;

    bamToFasta(subjectBam);
    
    return 0;   
}

