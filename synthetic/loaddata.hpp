#ifndef _LOADDATA_H
#define _LOADDATA_H

/**
 * Split a string into a vector of sub-strings, according to sperator specified by c
 */
std::vector<std::string> split(const std::string &str, const char *c);

/**
 * parse a line into key-value pairs
 */
bool parse_row(const std::string &line, std::tuple<unsigned int, unsigned int> &t);

/**
 * parse a file into a group of key-value pairs
 */
std::vector<std::tuple<unsigned int, unsigned int> > load_data(const char *file);

#endif
