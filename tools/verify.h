#ifndef __VERIFY_H__
#define __VERIFY_H_

class verify{
 public:
  verify(std::string, int32_t, std::string, std::vector<std::string>, std::string);
 

  
 private:
  const int32_t kmerSize_;
  const std::string sequence_;
  const std::string probandPath_;
  const std::string hashListPath_;
  const std::vector<std::string> controlPaths_;
  std::vector<std::string> sequenceKmers_;
  std::map<std::string, int32_t> probandKmers_;
  std::map<std::string, std::map<std::string, int32_t> > controlKmers_;

  void printKmerMap(const std::map<std::string, int32_t> &, const std::string &);
  const std::vector<std::string> kmerize();
  const std::map<std::string, int32_t> countKmers(const std::string &);
  const std::map<std::string, int32_t> getKmersFromHashList(const std::string &);
  const std::map<std::string, int32_t> filterHashListKmers(const std::map<std::string, int32_t> &);

};

#endif // __VERIFY_H__
