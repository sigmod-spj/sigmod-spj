#include <cstring>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>

using namespace std;


/**
 * Split a string into a vector of sub-strings, according to sperator specified by c
 * @param str input string
 * @param seperator list
 */
vector<string> split(const string &str, const char *c) {
  char *cstr, *p;
  vector<string> res;
  cstr = new char[str.size() + 1];
  strcpy(cstr, str.c_str());
  p = strtok(cstr, c);
  while(p != NULL) {
    res.emplace_back(p);
    p = strtok(NULL,c);
  }
  return res;
}

/**
 * parse a line into key-value tuples
 * @param line the input line to be parsed
 * @param t the tuple where the values are stored
 * @return true if the pair is generated, false if the line starts with # or does not have exactly 2 values
 */
bool parse_row(const string &line, tuple<unsigned int, unsigned int> &t) {
  if (line[0] == '#') {
    return false;
  } else {
    vector<string> tokens = split(line, " \t");
    if (tokens.size() != 2) {
      return false;
    } else {
      t = make_tuple(stoul(tokens[0]),stoul(tokens[1]));
      return true;
    }
  }
}

/**
 * parse a file into a group of key-value pairs
 * @param file filename
 * @return a vector of tuples representing the relation
 */
vector<tuple<unsigned int, unsigned int> > load_data(const char *file) {
  vector<tuple<unsigned int, unsigned int> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  while(getline(fin,line)){
    tuple<unsigned int, unsigned int>  t;
    if (parse_row(line, t)) {
      result.emplace_back(t);      
    }
  }
  return result;
}

