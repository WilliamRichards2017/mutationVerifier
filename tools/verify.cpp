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

//#include <cxxopts.hpp>
#include "../bin/externals/cxxopts/src/cxxopts_project/include/cxxopts.hpp"

struct region{
  std::string startChrom;
  int32_t startPos;
  std::string endChrom;
  int32_t endPos;
};

struct proband{
  std::vector<BamTools::BamAlignment> reads;
};

const std::map<std::string, int32_t> verify::getSequenceCountsFromKmerMap(const std::map<std::string, int32_t> & kmerMap){
  std::map<std::string, int32_t> sequenceMap;
  for (const auto & k : sequenceKmers_) {
    if(kmerMap.count(k) == 1){
      std::cout << "found kmer " << k << " with count " << kmerMap.at(k) << std::endl;
    }
    std::cout << k << " has count 0" << std::endl;
  }
  return sequenceMap;
}

void verify::printAllSequenceCounts(){
  std::map<std::string, int32_t> temp = verify::getSequenceCountsFromKmerMap(probandKmers_);
  for(const auto & kv : controlKmers_){
    std::map<std::string, int32_t> temp2 = verify::getSequenceCountsFromKmerMap(kv.second);
  }
}


verify::verify(std::string sequence, int32_t kmerSize, std::string probandPath, std::vector<std::string> controlPaths) : sequence_(sequence), kmerSize_(kmerSize), probandPath_(probandPath), controlPaths_(controlPaths){

  sequenceKmers_ = verify::kmerize();

  for(const auto & c : controlPaths_){
    std::cout << "control file is: " << c << std::endl;
    std::map<std::string, int32_t> controlKmer = verify::getKmersFromJhash(c);
    controlKmers_.insert({c, controlKmer});
  }

  probandKmers_ = verify::getKmersFromJhash(probandPath_);

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

const std::map<std::string, int32_t> verify::getKmersFromJhash(const std::string & jhashPath){
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


const std::vector<std::string> verify::kmerize(){
  int32_t kmercount = 0;

  std::vector<std::string> kmers;

  std::cout << "string to kmerize is: " << sequence_ << std::endl;
  
  while(kmercount + kmerSize_ <= sequence_.length()){
    std::string kmer = sequence_.substr(kmercount, kmerSize_);
    kmers.push_back(kmer);
    ++kmercount;

    
    std::cout << "pushing back kmer: " << kmer << std::endl;
    std::cout << "kmer count is: " << kmercount << std::endl;
    std::cout << "kmersize is: " << kmerSize_ << std::endl;
    std::cout << "sequence length is: " << sequence_.length() << std::endl;
    std::cout << "incremented pos is: " << kmercount+kmerSize_ << std::endl;
    
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
  v.printAllSequenceCounts();

  return 0;
}
