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

const std::map<std::string, int32_t> verify::getSequenceCountsFromKmerMap(const std::map<std::string, int32_t> & kmerMap){
  std::map<std::string, int32_t> sequenceMap;
  std::cout << "size of kmer map is: " << kmerMap.size() << std::endl;
  for (const auto & k : sequenceKmers_) {
    std::string revK = util::revComp(k);
    if(kmerMap.count(k) != 0){
      std::cout << "found kmer " << k << " with count " << kmerMap.at(k) << std::endl;
    }
    else if (kmerMap.count(revK) != 0){
      std::cout << "found kmer rev comp " << revK << " with count " << kmerMap.at(revK) << std::endl;
    }
    else{
      std::cout << k << " has count 0" << std::endl;
    }
  }
  return sequenceMap;
}

/*void verify::printAllSequenceCounts(){
  std::cout << "######################_PROBAND_###################################" << std::endl;
  std::map<std::string, int32_t> temp = verify::getSequenceCountsFromKmerMap(probandKmers_);
  std::cout << "##################################################################" << std::endl;
  for(const auto & kv : controlKmers_){
    std::map<std::string, int32_t> temp2 = verify::getSequenceCountsFromKmerMap(kv.second);
  }
  }*/

/*const std::string jhashToDumpPath(const std::string & jhashPath){
  std::string dumpPath = "";
  dumpPath += "/uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/resources/testData/";
  dumpPath += util::baseName(jhashPath);
  dumpPath += ".dump";
  return dumpPath;
  }*/

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

    //std::cout << "kmerCount[0], atoi(kmerCount[1]) is: " << kmerCount[0] << ", " << atoi(kmerCount[1].c_str()) << std::endl;
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
  //probandDumpPath_ = jhashToDumpPath(probandPath_);
  //std::cout << "proband dump path is " << probandDumpPath_ << std::endl;
  
  std::map<std::string, int32_t> probandKmerCounts = verify::countKmers(probandPath_);
  verify::printKmerMap(probandKmerCounts, probandPath_);
  for(const auto & c : controlPaths_){
    std::map<std::string, int32_t> controlKmerCounts = verify::countKmers(c);
    verify::printKmerMap(controlKmerCounts, c);
  }
  

  probandKmers_ = verify::getKmersFromJhash(probandDumpPath_);

}

const std::string & verify::getSequence(){
  return sequence_;
}

const std::string & verify::getProbandPath(){
  return probandPath_;
}

const std::vector<std::string> & verify::getControlPaths(){
  return controlPaths_;
}

const std::vector<std::string> & verify::getSequenceKmers(){
  return sequenceKmers_;
}

const std::map<std::string, std::map<std::string, int32_t> > & verify::getControlKmers(){
  return controlKmers_;
}


/*
  void bamToFasta(std::string bamFile){
  std::string command = "samtools bam2fq ";
  command += bamFile;
  command += " | seqtk seq -A - > temp/";
  command += util::baseName(bamFile);
  command += ".fa";

  std::cout << "running command " << command << std::endl;
  util::exec(command.c_str());
}
*/

/*void dumpJhash(std::string jhashFile){
  std::string jellyfishPath = "/uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";

  std::string histoCommand = "";
  histoCommand += jellyfishPath;
  histoCommand += " histo ";
  histoCommand += jhashFile;

  std::cout << "executing histo command " << histoCommand << std::endl;
  util::exec(histoCommand.c_str());

  std::string dumpCommand = "";
  dumpCommand += jellyfishPath;
  dumpCommand += " dump ";
  dumpCommand += jhashFile;
  dumpCommand += " > ";
  dumpCommand += "/uufs/chpc.utah.edu/common/home/u0401321/mutationVerifier/resources/testData/";
  dumpCommand += util::baseName(jhashFile);
  dumpCommand += ".dump";
 
  
  //  std::cout<< "executing dump command " << dumpCommand << std::endl;
  util::exec(dumpCommand.c_str());
  }*/

const std::map<std::string, int32_t> verify::getKmersFromJhash(const std::string & jhashPath){
  std::ifstream file(jhashPath);
  std::string line;

  std::cout << "reading in kmers from " << jhashPath << std::endl;

  std::map<std::string, int32_t> kmerMap;
  
  std::ifstream newfile(jhashPath);
  while(std::getline(newfile, line)){
    std::string countString = line.erase(0,1);
    int32_t count = atoi(countString.c_str());
    
    //std::cout << "count is: " << count << std::endl;
    std::getline(newfile, line);
    std::string kmer = line;

    //std::cout << "inserting element into map: " << kmer << ":" << count << std::endl;
    kmerMap.insert({kmer, count});
    //std::cout << "kmer is " << kmer << std::endl;
  }
  return kmerMap;
}


const std::vector<std::string> verify::kmerize(){
  int32_t kmercount = 0;

  std::vector<std::string> kmers;

  std::cout << "string to kmerize is: " << sequence_ << std::endl;
  
  while(kmercount + kmerSize_ <= sequence_.length()){
    std::string kmer = sequence_.substr(kmercount, kmerSize_);
    kmers.push_back(kmer);
    ++kmercount;

    /*
    std::cout << "pushing back kmer: " << kmer << std::endl;
    std::cout << "kmer count is: " << kmercount << std::endl;
    std::cout << "kmersize is: " << kmerSize_ << std::endl;
    std::cout << "sequence length is: " << sequence_.length() << std::endl;
    std::cout << "incremented pos is: " << kmercount+kmerSize_ << std::endl;
    */
    
  }
  return kmers;
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
    //dumpJhash(probandPath);
  }
  else{
    std::cout << "Please provide a proband Jhash file with [-p|--proband]" << std::endl;
    std::cout << "Exiting run with non-zero exit status" << std::endl;
    exit (EXIT_FAILURE);
  }
  
  std::vector<std::string> controlPaths;
  if(result.count("c")){
    controlPaths = result["c"].as<std::vector<std::string>>();
    for(const auto & c : controlPaths){
      //dumpJhash(c);
    }
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
  //v.printAllSequenceCounts();

  return 0;
}
