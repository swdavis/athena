// spectrum_reader.cpp

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "spectrum_reader.hpp"
#include "../athena.hpp"

void ReadSpectrumToCDF(const std::string& filename, std::vector<Real>& wl, std::vector<Real>& cdf, Real& itot) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("ReadSpectrumToCDF: cannot open file " + filename);
  }

  // clear any values that might already be in wl and cdf
  wl.clear();
  cdf.clear();

  std::vector<Real> ilam;
  std::string line;
  int line_number = 0;
  while (std::getline(file, line)) {
    ++line_number;

    // skip empty lines and comments
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    Real wavelength, intensity;

    // check for lines that don't read nicely into two doubles
    if (!(iss >> wavelength >> intensity)) {
      throw std::runtime_error("ReadSpectrumToCDF: malformed line " + std::to_string(line_number) + " in input file '" + filename + "'");
    }

    wl.push_back(wavelength);
    ilam.push_back(intensity);

  }
  file.close();

  const int nrows = wl.size();
  if (nrows == 0) {
    throw std::runtime_error("ReadSpectrumToCDF: no valid data in file " + filename);
  }
  if (nrows < 2) {
    throw std::runtime_error("ReadSpectrumToCDF: need at least 2 data points for interpolation, found " + std::to_string(nrows));
  }

  // compute CDF as int dlambda*ilam
  // assumes wavelengths are uniformly spaced and increasing
  Real dlambda = wl[1] - wl[0];
  //printf("dlambda=%g\n", dlambda);
  cdf.push_back(0.);
  for (int i = 1; i < nrows; ++i) {
    cdf.push_back(cdf[i-1] + ilam[i]*dlambda);
  }
  
  // normalize: CDF runs from 0 to 1
  itot = cdf[nrows-1];
  printf("ReadSpectrumToCDF: itot = %g [erg/cm^2/s]\n", itot);
  if (itot == 0.) {
    throw std::runtime_error("ReadSpectrumToCDF: cdf norm is zero");
  }
  for (int i=0; i<nrows; ++i) {
    cdf[i] = cdf[i] / itot;
    //printf("i=%d, wl=%g, ilam=%g, cdf=%g\n", i, wl[i], ilam[i], cdf[i]);
  }

} // end ReadSpectrumToCDF
