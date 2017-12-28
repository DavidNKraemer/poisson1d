#ifndef _UTILS_H
#define _UTILS_H

#include <fstream>
#include <exception>
#include <vector>

namespace ams562 {
  void write3v2txt(std::string &filename, std::vector<double> xdata,
    std::vector<double> ydata, std::vector<double> ydata_exact);
}

#endif
