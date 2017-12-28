#include "utils.h"

namespace ams562 {
  void write3v2txt(std::string &filename, std::vector<double> xdata,
      std::vector<double> ydata, std::vector<double> ydata_exact) {

    std::ofstream outfile(filename);

    if (outfile.is_open()) {

      for (unsigned int i = 0; i < xdata.size(); i++){
        outfile << xdata[i] << "," << ydata[i] << "," << ydata_exact[i] << std::endl;
      }

      outfile.close();
    }

    else {
      std::string err { "[Error] Cannot open: " + filename + " to write!"};
      throw err;
    }
  }
}
