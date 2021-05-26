#include <fstream>
#include <iostream>
#include <omp.h>

#include "lib/SboxTools.hpp"

int main() {
  const size_t outdim = 1;
  const size_t dim = 5;
  std::string filename = "../SboxData/F5b.txt";
  std::ifstream file(filename);
  auto functions=
      SboxTools::ReadFile(file, dim, outdim);
  SboxTools::csv_line::PrintHeader(std::cout);
  for (const auto &x_raw: functions) {
    VectorialBooleanFunction x(x_raw, dim, outdim);
    auto line = SboxTools::Analyse<SboxTools::csv_line>(x);
    {
      std::cout << line;
    }
  }
}
