#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdint.h>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

#include <regex.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "util.h"

//#include <cxxopts.hpp>
#include "/uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/bin/externals/cxxopts/src/cxxopts_project/include/cxxopts.hpp"

struct region{
  std::string startChrom;
  int32_t startPos;
  std::string endChrom;
  int32_t endPos;
};

struct proband{
  std::vector<BamTools::BamAlignment> reads;
};

void bamToFasta(std::string bamFile){
  std::string command = "samtools bam2fq ";
  command += bamFile;
  command += " | seqtk seq -A - > temp/";
  command += util::baseName(bamFile);
  command += ".fa";

  std::cout << "running command " << command << std::endl;
  util::exec(command.c_str());
}

void runJelly(std::string fastaFile, int32_t kmerSize, int32_t threads){
  std::string jellyfishPath = "/uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";
  std::string command = "";
  command += jellyfishPath;
  command += " count -m ";
  command += std::to_string(kmerSize);
  command += " -s 100M -t ";
  command += std::to_string(threads);
  command += " -o /uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/temp/fasta.jf";
  command += " -C ";
  command += fastaFile;

  std::cout << "executing command " << command << std::endl;
  util::exec(command.c_str());

  std::string histoCommand = "";
  histoCommand += jellyfishPath;
  histoCommand += " histo /uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/temp/fasta.jf";

  std::cout << "executing histo command " << histoCommand << std::endl;
  util::exec(histoCommand.c_str());
  
  std::string dumpCommand = "";
  dumpCommand += jellyfishPath;
  dumpCommand += " dump /uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/temp/fasta.jf > /uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/temp/fasta.Jhash";
  
  std::cout << "executing dump command " << dumpCommand << std::endl;
  util::exec(dumpCommand.c_str());
  
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

  proband p;

  while(reader.GetNextAlignment(al)){
    p.reads.push_back(al);
    std::cout << "found read in region" << std::endl;
  }

}

std::map<std::string, int32_t> jhashToMap(std::string jhashPath){
  std::ifstream file(jhashPath);
  std::string line;

  std::map<std::string, int32_t> kmerMap;
  
  std::ifstream newfile(jhashPath);
  while(std::getline(newfile, line)){
    std::string countString = line.erase(0,1);
    int32_t count = atoi(countString.c_str());
    
    std::cout << "count is: " << count << std::endl;
    std::getline(newfile, line);
    std::string kmer = line;

    //std::cout << "inserting element into map: " << kmer << ":" << count << std::endl;
    kmerMap.insert({kmer, count});
    std::cout << "kmer is " << kmer << std::endl;
  }
  return kmerMap;
}


std::map<std::string, int32_t> getUniqueHashes(std::string probandUniqueHashList, std::string regionJhash){
  //void getUniqueHashes(std::string probandUniqueHashList, std::string regionJhash){
  return jhashToMap(regionJhash);
}

int main(int argc, char* argv[]){
  cxxopts::Options options(argv[0], "Verifies if proband dna is unique to proband sample");

  options.add_options()
    ("s,sequence", "DNA sequence to verify uniqueness", cxxopts::value<std::string>())
    ("p,probandJhash", "proband Jhash file", cxxopts::value<std::string>())
    ("help", "Print help")
    ("c,controlJhashes", "comma-seperated list of control Jhash files", cxxopts::value<std::string>());
  
  auto result = options.parse(argc, argv);

  if(result.count("sequence")){
    std::cout << "sequence is: " << result["s"].as<std::string>() << std::endl;
    }
  
  std::string seq = result["s"].as<std::string>();
  std::cout << "parsed sequence is: " << seq << std::endl;
  
  return 0;
}
