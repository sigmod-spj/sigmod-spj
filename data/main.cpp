/**
 * Data Generator for synthetic datasets
 */
#include <fstream>
#include <iostream>
#include <string>  /* getline */
#include <math.h>  /* pow */
using namespace std;

/**
 * Generates uniform data
 * Tuples are in the form of "i j" where both i and j ranges from 1 to 1000 by default
 */
void generate_uniform_data(unsigned int D = 1000, unsigned int freq = 1000) {
  char filename[] = "uniform.dat";
  ofstream outfile;
  outfile.open(filename);
  unsigned int lines = 0;
  for (unsigned int i = 1; i <= D; ++i) {
    for (unsigned int j = 1; j <= freq; ++j) {
      outfile << i << ' ' << j << endl;
      lines += 1;
    }
  }
  outfile.close();
  cout << lines << " lines of uniform data written to \'" << filename << "\'." << endl;
}

/**
 * Generates Zipf data
 * Tuples are in the form of "i j" where both i ranges from 1 to D and j ranges from 1 to freq_i
 */
void generate_zipf_data(unsigned int D = 50000, float z = 1.1) {
  char filename[] = "zipf.dat";
  ofstream outfile;
  outfile.open(filename);
  unsigned int lines = 0;

  for (unsigned int i = 1; i <= D; ++i) {
    unsigned int freq = round( pow( static_cast<double>(D)/i , z ) );
    for (unsigned int j = 1; j <= freq; ++j) {
      outfile << i << ' ' << j << endl;
      lines += 1;
    }
  }
  outfile.close();
  cout << lines << " lines of zipfian data written to \'" << filename << "\'." << endl;
}

int main() {
  string input;
  cout << "Input the type of data to generate [default = all]:" << endl;
  cout << "1. Uniform Data" << endl;
  cout << "2. Zipfian Data" << endl;
  cout << "Type: ";
  int c = 0;
  getline(cin,input);
  if (!input.empty()) {
    c = stoi(input);
  }

  if ((c == 0) || (c == 1)) {
    cout << "-- Generating uniform data --" << endl;
    cout << "Input D [default = 1000]: ";
    unsigned int D = 1000;
    getline(cin,input);
    if (!input.empty()) {
      D = stoi(input);
    }

    cout << "Input uniform frequency [default = 1000]: ";
    unsigned int freq = 1000;
    getline(cin,input);
    if (!input.empty()) {
      freq = stoi(input);
    }

    generate_uniform_data(D, freq);

  }

  if ((c == 0) || (c == 2)) {
    cout << "-- Generating zipf data --" << endl;
    cout << "Input D [default = 50000]: ";
    unsigned int D = 50000;
    getline(cin,input);
    if (!input.empty()) {
      D = stoi(input);
    }

    cout << "Input Skewness Parameter [default = 1.1]: ";
    float z = 1.1;

    getline(cin,input);
    if (!input.empty()) {
      z = stof(input);
    }

    generate_zipf_data(D, z);

  }

}
