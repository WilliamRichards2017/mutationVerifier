#ifndef __VERIFY_H__
#define __VERIFY_H_

class verify{
 public:
  verify(std::string, int32_t, std::string, std::vector<std::string>);
 
  const std::string & getSequence();
  const std::string & getProbandPath();
  const std::vector<std::string> & getControlPaths();
  const std::vector<std::string> & getSequenceKmers();
  const std::map<std::string, std::map<std::string, int32_t> > & getControlKmers();
  
 private:
  const int32_t kmerSize_;
  const std::string sequence_;
  const std::string probandPath_;
  const std::vector<std::string> controlPaths_;
  const std::vector<std::string> sequenceKmers_;
  std::map<std::string, std::map<std::string, int32_t> > controlKmers_;

  const std::vector<std::string> kmerize();
  const std::map<std::string, int32_t> getKmersFromJhash(const std::string &);
};

#endif // __VERIFY_H__