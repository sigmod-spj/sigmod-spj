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

std::vector<std::tuple<unsigned int, unsigned int, int> > load_graph(const char *file);

bool parse_store_sales(const std::string &line, std::tuple<unsigned int, unsigned int, double, double> &t);

std::vector<std::tuple<unsigned int, unsigned int, double, double> > load_store_sales(const char *file);

std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > load_date(const char *file);

std::vector<std::tuple<unsigned int, std::string, std::string> > load_item(const char *file);

std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > load_sales(const char *file, bool is_catelog);

std::vector<std::tuple<unsigned int, unsigned int, unsigned short> > load_principals(const char *file);

std::vector<std::tuple<unsigned int, double> > load_ratings(const char *file);

#endif
