#ifndef CSVWRITER_HPP
#define CSVWRITER_HPP
#include <fstream>
#include <sstream>
#include <lsSmartPointer.hpp>
// A simple CSV writer class
template <class NumericType> class CSVWriter {
  std::ofstream file;

  template <class Iterator>
  static std::string join(Iterator begin, Iterator end,
                          const std::string &separator = ",") {
    std::ostringstream ostr;
    if (begin != end)
      ostr << *begin++;
    while (begin != end)
      ostr << separator << *begin++;
    return ostr.str();
  }

public:
  CSVWriter(std::string passedFilename) : file(passedFilename) {}

  void writeRow(const std::vector<NumericType> &data) {
    file << join(data.cbegin(), data.cend()) << "\n";
  }

  void writeRow(lsSmartPointer<std::vector<NumericType>> data) {
    if (data == nullptr)
      return;
    file << join(data->cbegin(), data->cend()) << "\n";
  }

  void writeLine(const std::string &line) { file << line << "\n"; }

  ~CSVWriter() { file.close(); }
};
#endif