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

//TODO: implement indexing bam
void indexBam(std::string bamFile){
}

void writeBamRegion(region reg, std::string bamPath){
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(bamPath)){
    std::cout << "could not open input Bam Path in writeBamRegion for " << bamPath << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }
  
  reader.LocateIndex();
  
  if(!reader.HasIndex()){
    std::cout << "Index for " << bamPath << "could not be opened in writeBamRegion" << std::endl;
    std::cout << "Exiting run with non-zero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }
  
  int32_t startRefID = reader.GetReferenceID(reg.startChrom);
  int32_t endRefID = reader.GetReferenceID(reg.endChrom);

  std::cout << "trying to set region for coords " << startRefID << ", " << reg.startPos << ", " << endRefID << ", " << reg.endPos << std::endl;

  BamTools::BamRegion region = {startRefID, reg.startPos, endRefID, reg.endPos};

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " << startRefID << ", " << reg.startPos << ", " << endRefID << ", " << reg.endPos << std::endl;
  }
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
    writeBamRegion(reg, subjectBam);
    
    return 0;   
}

