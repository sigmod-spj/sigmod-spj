#include <cstring>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>

using namespace std;


/**
 * split a string into a vector of sub-strings, according to sperator specified by c
 * empty-column aware
 */

vector<string> split(const string& src, const char& delim) {
  vector<string> vec;
  int src_len = src.length();
  int find_cursor = 0;
  int read_cursor = 0;

  if (src_len <= 0) return vec;

  while (read_cursor < src_len){

    find_cursor = src.find(delim, find_cursor);

    if (-1 == find_cursor) {
      if (read_cursor <= 0) {
	return vec;
      }

      if (read_cursor < src_len){
	vec.push_back(src.substr(read_cursor, src_len - read_cursor));
	return vec;
      }
    } else if (find_cursor == read_cursor) {
      vec.push_back(string(""));
    } else {
      vec.push_back(src.substr(read_cursor, find_cursor - read_cursor));
    }
    
    read_cursor = ++find_cursor;
    if (read_cursor == src_len){
      vec.push_back(string(""));
    } 
  }//end while()

  return vec;
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
    vector<string> tokens = split(line, '|');
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

/**
 * parse a file into a group of key-value pairs
 * @param file filename
 * @return a vector of tuples representing the relation
 */
vector<tuple<unsigned int, unsigned int, int> > load_graph(const char *file) {
  vector<tuple<unsigned int, unsigned int, int> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  while(getline(fin,line)){
    vector<string> tokens = split(line, ',');
    unsigned int source = (tokens[0].empty())? 0 : stoul(tokens[0]);
    unsigned int target = (tokens[1].empty())? 0 : stoul(tokens[1]);
    int rating = (tokens[2].empty())? 0 : stoi(tokens[2]);
    //    double time = (tokens[3].empty())? 0 : stod(tokens[3]);
    result.emplace_back
      (make_tuple(source, target, rating));      
  }
  return result;
}


/**
 * load the store_sales table
 * @param file filename
 * @return a vector of tuples representing the relation
 * columns: 100 * ss_list_price, ss_quantity, ss_coupon_amt, ss_wholesale_cost
 */
vector<tuple<unsigned int, unsigned int, double, double> > load_store_sales(const char *file) {
  vector<tuple<unsigned int, unsigned int, double, double> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  while(getline(fin,line)){
    vector<string> tokens = split(line, '|');
    double ss_list_price = (tokens[12].empty())? 0 : stod(tokens[12]);
    unsigned int ss_quantity = (tokens[10].empty())? 0 : stoul(tokens[10]);
    double ss_coupon_amt = (tokens[19].empty())? 0 : stod(tokens[19]);
    double ss_wholesale_cost = (tokens[11].empty())? 0 : stod(tokens[11]);
    result.emplace_back
        (make_tuple(static_cast<unsigned int>(100 * ss_list_price),
		 ss_quantity,
		 ss_coupon_amt,
		 ss_wholesale_cost));
  }
  return result;
}

/**
 * load the date table
 * @param file filename
 * @return a vector of tuples representing the relation
 * columns: d_date_sk, d_moy, d_year
 */
vector<tuple<unsigned int, unsigned int, unsigned int> > load_date(const char *file) {
  vector<tuple<unsigned int, unsigned int, unsigned int> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  while(getline(fin,line)){
    vector<string> tokens = split(line, '|');
    unsigned int d_date_sk = stoul(tokens[0]);
    unsigned int d_moy = stoul(tokens[8]);
    unsigned int d_year = stoul(tokens[6]);

    result.emplace_back
        (make_tuple(d_date_sk,
                    d_moy,
                    d_year));
  }
  return result;
}

/**
 * load the item table
 * @param file filename
 * @return a vector of tuples representing the relation
 * columns: i_item_sk, i_category_id, i_category, i_class_id, i_class
 */
vector<tuple<unsigned int, string, string> > load_item(const char *file) {
  vector<tuple<unsigned int, string, string> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  while(getline(fin,line)){
    vector<string> tokens = split(line, '|');

    unsigned int i_item_sk = stoul(tokens[0]);
    string i_category = tokens[12];
    string i_class = tokens[10];

    result.emplace_back
        (make_tuple(i_item_sk,
                    i_category,
                    i_class));
  }
  return result;
}

/**
 * load the sales table
 * @param file filename
 * @return a vector of tuples representing the relation
 * columns: customer_sk, sold_date_sk, item_sk
 */
vector<tuple<unsigned int, unsigned int, unsigned int> > load_sales(const char *file, bool is_catelog) {
  vector<tuple<unsigned int, unsigned int, unsigned int> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  while(getline(fin,line)){
    vector<string> tokens = split(line, '|');

    unsigned int customer_sk;
    unsigned int sold_date_sk;
    unsigned int item_sk;

    if (is_catelog) {
      customer_sk = (tokens[3].empty())?(-1):stoul(tokens[3]);
      sold_date_sk = (tokens[0].empty())?(-1):stoul(tokens[0]);
      item_sk = stoul(tokens[15]);
    } else {
      customer_sk = (tokens[4].empty())?(-1):stoul(tokens[4]);
      sold_date_sk = (tokens[0].empty())?(-1):stoul(tokens[0]);
      item_sk = stoul(tokens[3]);
    }


    result.emplace_back
        (make_tuple(customer_sk,
                    sold_date_sk,
                    item_sk));
  }
  return result;
}

/**
 * @param file filename
 * columns: nconst, tconst, category
 */
vector<tuple<unsigned int, unsigned int, unsigned short> > load_principals(const char *file) {
  vector<tuple<unsigned int, unsigned int, unsigned short> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  getline(fin,line); // skip the first line
  while(getline(fin,line)){
    vector<string> tokens = split(line, '\t');

    unsigned int tconst = stoul(tokens[0].erase(0,2)); // drop frist two 
    unsigned int nconst = stoul(tokens[2].erase(0,2));

    
    unsigned short category = 0;
    if (tokens[3].compare("director") == 0) {
      category = 1;
    } else if (tokens[3].compare("actor") == 0) {
      category = 2;
    } else if (tokens[3].compare("actress") == 0) {
      category = 3;
    } else if (tokens[3].compare("self") == 0) {
      category = 4;
    } else if (tokens[3].compare("cinematographer") == 0) {
      category = 5;
    } else if (tokens[3].compare("composer") == 0) {
      category = 6;
    } else if (tokens[3].compare("editor") == 0) {
      category = 7;
    } else if (tokens[3].compare("producer") == 0) {
      category = 8;
    } else if (tokens[3].compare("writer") == 0) {
      category = 9;
    } 
    result.emplace_back
      (make_tuple(nconst, tconst, category));
  }
  return result;
}

/**
 * @param file filename
 * columns: nconst, tconst, category
 */
vector<tuple<unsigned int, double> > load_ratings(const char *file) {
  vector<tuple<unsigned int, double> > result;
  ifstream fin;
  fin.open(file, ifstream::in);
  if(!fin.good())
    return result;
  string line;
  getline(fin,line); // skip the first line
  while(getline(fin,line)){
    vector<string> tokens = split(line, '\t');

    unsigned int tconst = stoul(tokens[0].erase(0,2)); // drop frist two 
    double rating = stod(tokens[1]);

    result.emplace_back
      (make_tuple(tconst, rating));
  }
  return result;
}

