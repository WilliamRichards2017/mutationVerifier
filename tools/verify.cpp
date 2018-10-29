#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include <regex.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "util.h"
#include "verify.h"

#include "../bin/externals/cxxopts/src/cxxopts_project/include/cxxopts.hpp"

const std::vector<std::string> verify::kmerize(){
  int32_t kmercount = 0;

  std::vector<std::string> kmers;

  std::cout << "string to kmerize is: " << sequence_ << std::endl;

  while(kmercount + kmerSize_ <= sequence_.length()){
    std::string kmer = sequence_.substr(kmercount, kmerSize_);
    kmers.push_back(kmer);
    ++kmercount;
  }

  return kmers;
}

const std::map<std::string, int32_t> verify::countKmers(const std::string & jhashPath){
  std::string jellyfishPath = "/uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";
  
  std::map<std::string, int32_t> ret;
  for (const auto & kmer : sequenceKmers_){

    std::string cmd = jellyfishPath + " query " + jhashPath + " " + kmer;
    //std::cout << "executing command: " << cmd << std::endl;
    
    std::string queryOutput = util::exec(cmd.c_str());
    //std::cout << "command output is: " << queryOutput << std::endl;
    std::istringstream iss(queryOutput);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)),
				       std::istream_iterator<std::string>());

    if(kmerCount.size() == 2){
      ret.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
      }
  }
  return ret;
}

void verify::printKmerMap(const std::map<std::string, int32_t> & kmerMap, const std::string & originPath ){
  std::cout << "Looking for mutations in Jhash file: " << originPath << std::endl;
  for(const auto & kv : kmerMap){
    if(kv.second < 2){
      std::cout << "found mutation kmer: " << kv.first << "  with count: " << kv.second << std::endl;
    }
  }
}

verify::verify(std::string sequence, int32_t kmerSize, std::string probandPath, std::vector<std::string> controlPaths) : sequence_(sequence), kmerSize_(kmerSize), probandPath_(probandPath), controlPaths_(controlPaths){

  sequenceKmers_ = verify::kmerize();
  std::map<std::string, int32_t> probandKmerCounts = verify::countKmers(probandPath_);
  verify::printKmerMap(probandKmerCounts, probandPath_);
  for(const auto & c : controlPaths_){
    std::map<std::string, int32_t> controlKmerCounts = verify::countKmers(c);
    verify::printKmerMap(controlKmerCounts, c);
  }
}

int main(int argc, char* argv[]){
  cxxopts::Options options(argv[0], "Verifies if proband dna is unique to proband sample");

  options.add_options()
    ("help", "Print help")
    ("s,sequence", "DNA sequence to verify uniqueness", cxxopts::value<std::string>())
    ("p,probandJhash", "proband Jhash file", cxxopts::value<std::string>())
    ("l,length", "length of kmer", cxxopts::value<int32_t>())
    ("c,controlJhashes", "comma-seperated list of control Jhash files", cxxopts::value<std::vector<std::string>>());
  
  auto result = options.parse(argc, argv);


  std::string sequence;
  if(result.count("sequence")){
    sequence = result["s"].as<std::string>();
    std::cout << "sequence is: " << sequence << std::endl;
  }
  else{
    std::cout << "Please provide a sequence with [-s|--sequence]" << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit (EXIT_FAILURE);
  }
  std::string probandPath;
  if(result.count("p")){
    probandPath = result["p"].as<std::string>();
    std::cout << "proband file is: " << probandPath << std::endl;
  }
  else{
    std::cout << "Please provide a proband Jhash file with [-p|--proband]" << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit (EXIT_FAILURE);
  }
  std::vector<std::string> controlPaths;
  if(result.count("c")){
    controlPaths = result["c"].as<std::vector<std::string>>();
  }
  else{
    std::cout << "Please provide atleast one control Jhash file  with [-c|--control]" << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit (EXIT_FAILURE);
  }
  

  int32_t kmerSize;
  if(result.count("l")){
    kmerSize = result["l"].as<int32_t>();
    std::cout << "kmer length is: " << kmerSize << std::endl;
  }
  else{
    std::cout << "Please provide a kmer length  with [-l|--length]" << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit (EXIT_FAILURE);
  }

  verify v = {sequence, kmerSize, probandPath, controlPaths};
 
  return 0;
}
