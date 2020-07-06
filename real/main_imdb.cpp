#include <iostream>     // std::cout
#include <fstream>
#include <algorithm>    // std::shuffle
#include <iomanip>      // std::setprecision
#include <vector>
#include <tuple>
#include <unordered_map>
#include <map>
#include <cstdlib>      // std::RANDMAX
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <limits>
#include "hasher.hpp"
#include "loaddata.hpp"

using namespace std;
using namespace std::chrono;

typedef vector<tuple<unsigned int, unsigned int, unsigned short> > TABLE_TYPE; 
typedef vector<tuple<unsigned int, unsigned int, unsigned short, double> > BINARY_JOIN_TYPE;
typedef vector<tuple<unsigned int, unsigned int, unsigned short, double> > RELATION_TYPE; 

// Filter actors(2) and actress(3) and rating >= param
template <class K, class... V>
bool tuple_filter(const tuple<K, V...>& t, int param) {
  return (get<2>(t) == 2 || get<2>(t) == 3) && 
    get<3>(t) >= param;
}

struct sketch_unit {
  unsigned int max_freq;
  unsigned int current_freq;
  bool was_complete;

  vector<double> prob;
  vector<double> temp_prob;

  TABLE_TYPE table;
  BINARY_JOIN_TYPE binary_join;
  
  sketch_unit(): max_freq(0), current_freq(0), was_complete(false) {}
  ~sketch_unit() {
    prob.clear();
    temp_prob.clear();
    table.clear();
    binary_join.clear();
  }

  double estimate(int param) const {
    bool success = false;
    double p = 1.0;
    for (unsigned int i = 0; i < binary_join.size(); ++i) {
      if (tuple_filter(binary_join[i], param)) {
  	if (!was_complete || success) { // light case or 2 successes or ss too small
  	  return 1.0;
  	} else {
  	  success = true;
  	  p = prob[i];
  	}
      }
    }
    if (success) { // 1 success
      return 1.0 / sqrt(p);
    } else { // 0 success
      return 0.0;
    }
  }
  
  void commit_table() {
    current_freq = 0;
    max_freq = max_freq * 10;
  }
  void commit_binary_join() {
    if (!was_complete && current_freq > max_freq) {
      was_complete = true;
      prob = vector<double>(binary_join.size(), max_freq / static_cast<double>(current_freq));
    }
    table.clear();
    current_freq = 0;
  }
  
  void add_table(const tuple<unsigned int, unsigned int, unsigned short> &s) {
    table.emplace_back(s);
    current_freq += 1;
    max_freq += 1;
  }

  // abcxy
  void add_binary_join(const tuple<unsigned int, unsigned int, unsigned short> &s,
		       const tuple<unsigned int, double> &d) {
    if (current_freq < max_freq) {
      binary_join.emplace_back
	(make_tuple(get<0>(s), get<1>(s), get<2>(s), get<1>(d)));
    } else {
      unsigned int idx = rand() % current_freq;
      if (idx < max_freq) {
    	binary_join[idx] = make_tuple(get<0>(s), get<1>(s), get<2>(s), get<1>(d));
      }
    }
    current_freq += 1;
  }
  
};

ostream& operator<<(ostream &os, sketch_unit s) {
  os << "MAX FREQUENCY: " << s.max_freq << endl
     << "CURRENT FREQUENCY: " << s.current_freq << endl
     << "WAS COMPLETE: " << s.was_complete << endl
     << "TABLE SIZE: " << s.table.size() << endl
     << "BINARY JOIN SIZE: " << s.binary_join.size() << endl
     << "PROB: ";
  for (const auto& p : s.prob) {
    os << p << ' ';
  }
  os << endl;
  return os;
}


template <class K, class... V>
vector<tuple<K, V...> > filter_table(const vector<tuple<K, V...> >& table, int i) {
  RELATION_TYPE result;
  for (const auto& t : table) {
    if (tuple_filter(t, i)) {
      result.emplace_back(t);
    }
  }
  return result;
}

// Warning: get_distinct_count will sort the data in place
template <class K, class... V>
unsigned int get_distinct_count(vector<tuple<K, V...> > &data) {
  if (data.empty()) {
    return 0;
  }
  sort(data.begin(), data.end());
  unsigned int d = 1;
  K current(get<0>(*data.cbegin()));
  for (auto it = ++data.cbegin(); it != data.cend(); ++it) {
    if (current != get<0>(*it)) {
      d += 1;
      current = get<0>(*it);
    }
  }
  return d;
}

/**
 * Solver using algorithm 1
 * @param data the input table
 * @param D the NDV of the key of data, auto-converted to double to avoid overflow
 * @param n sample size
 * @param M the M returned by algo. 1
 * @param p_i the parameters returnd by algo. 1 (We set p_i to 0 if tau_i==0)
 * @return the duration of the algorithm
 */
template <class KK, class... V>
duration<double> solve
(vector<tuple<KK, V...> > &data,
 double D,
 unsigned int n,
 unsigned int &M,
 unordered_map<KK, double> &p_i) {
  if (data.empty()) {
    return duration<double>();
  }
  sort(data.begin(), data.end());

  // Collect Degree Information
  vector<pair<unsigned int, KK> > N_i;
  unsigned int N_a = 1;
  KK current(get<0>(*data.cbegin()));
  for (auto it = ++data.cbegin(); it != data.cend(); ++it) {
    if (current == get<0>(*it)) {
      N_a += 1;
    } else {
      N_i.emplace_back(make_pair(N_a, current));
      current = get<0>(*it);
      N_a = 1;
    }
  }
  N_i.emplace_back(make_pair(N_a, current));

  // Sort value by degrees
  sort(N_i.begin(), N_i.end()); // increasing

  // find param starts
  high_resolution_clock::time_point start = high_resolution_clock::now();

  double min_error = std::numeric_limits<double>::max();
  double min_M = 0;
  double min_K = 0;

  double error = D * D;
  unsigned int M_0 = 0;

  double sum = 0.0;

  // Find M_0
  for (auto i = 1; i <= D; ++i) {
    sum += N_i[i-1].first;
    if (n >= sum) {
      M_0 = i;
    } else {
      sum -= N_i[i-1].first;
      break;
    }
  } // Now M_0 is the largest integer that n >= sum_1^M0 N_i

  unsigned int K = M_0; // sum is sum_1^K N_i
  double N_k = N_i[K-1].first;
  double sum_sqrt = -sqrt(N_k); // sum_sqrt will be sum_K+1^M sqrt(N_i)

  for (M = M_0; M <= D; ++M) {
    sum_sqrt += sqrt(N_i[M-1].first);
    while (K > 0 && n - sum <= sqrt(N_k) * sum_sqrt) { // K is too large
      sum -= N_k;
      sum_sqrt += sqrt(N_k);
      K = K-1;
      if (K > 0) {
        N_k = N_i[K-1].first;
      } else {
        N_k = 0;
      }
    }

    error = (D-M)*(D-M) + sum_sqrt * sum_sqrt / (n - sum) - M + K;

    if (error < min_error) {
      min_error = error;
      min_M = M;
      min_K = K;
    }
  }
  M = min_M;
  K = min_K;
  sum = 0.0;
  sum_sqrt = 0.0;

  for (int i = 0; i < K; ++i) {
    sum += N_i[i].first;
    p_i.emplace(N_i[i].second, 1.0);
  } // 1...K

  for (int i = K; i < M; ++i) {
    sum_sqrt += sqrt(N_i[i].first);
  }
  for (int i = K; i < M; ++i) {
    p_i.emplace(N_i[i].second, (n - sum) / sum_sqrt / sqrt(N_i[i].first));
  } // K+1 ... M

  for (int i = M; i < D; ++i) {
    p_i.emplace(N_i[i].second, 0.0);
  }

  cerr << "Parameter K = " << K << ", M = " << M << endl;
  high_resolution_clock::time_point end = high_resolution_clock::now();
  return duration_cast<duration<double>>(end - start);
}

template <class K, class... V>
vector<tuple<K, V...> > collect_uds_sketch
(const vector<tuple<K, V...> > &data,
 unsigned int n, unsigned int tau, unsigned int &l) {
  unordered_map<K, vector<tuple<K, V...> > > S;
  unsigned int S_size = 0;
  unordered_map<K, unsigned int> c;
  die_hasher h;
  l = 0;
  for (const auto& t : data) {
    K v = get<0>(t);
    unsigned int die_level = h(v);
    if (die_level >= l) {
      if (S.find(v)==S.end()) {
	S.emplace(v, RELATION_TYPE());
	S[v].emplace_back(t);
	c[v] = 1;
	S_size += 1;
      } else if (c[v] <= tau-1) {
	S[v].emplace_back(t);
	c[v] += 1;
	S_size += 1;
      } else {
	unsigned int idx = rand() % c[v];
	if (idx < tau) {
	  S[v][idx] = t;
	}
	c[v] += 1;
      }

      if (S_size >= n) {
	for (auto map_it = S.begin(); map_it != S.end();) {
	  K v = map_it -> first;
	  if (h(v) == l) {
	    S_size -= (map_it -> second).size();
	    map_it = S.erase(map_it);
	    c.erase(v);
	  } else {
	    ++map_it;
	  }
	}
	++l;
      }
    }
  }
  vector<tuple<K, V...> > result;
  for (const auto& map_t : S) {
    result.insert(result.end(),map_t.second.begin(),map_t.second.end());
  }
  return result;
}

template <class K, class... V>
unordered_map<K, vector<tuple<K, V...> > > collect_wds_sketch
(const vector<tuple<K, V...> > &data,
 const unordered_map<K, double>& p_i) {
  unordered_map<K, vector<tuple<K, V...> > > S;
  unsigned int S_size = 0;
  unordered_map<K, unsigned int> c;
  double_hasher h;

  for (const auto& t : data) {
    K a = get<0>(t);
    double v = h(a);
    if (v <= p_i.find(a) -> second) {
      if (S.find(a)!=S.end()) {
	S[a].emplace_back(t);
	c[a] += 1;
	S_size += 1;
      } else {
	S.emplace(a, RELATION_TYPE());
	S[a].emplace_back(t);
	c[a] = 1;
	S_size += 1;
      }
    }
  }

  return S;
}

template <class K, class... V>
bool tuple_less_key(const tuple<K, V...> &x, const K &y) {
  return get<0>(x) < y;
}

template <class K, class... V>
bool key_less_tuple(const K &x, const tuple<K, V...> &y) {
  return x < get<0>(y);
}

inline void set_timer(std::chrono::high_resolution_clock::time_point &timer) {
  timer = std::chrono::high_resolution_clock::now();
}

inline std::chrono::duration<double> get_timer(const std::chrono::high_resolution_clock::time_point &timer) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - timer);
}

int main(int argc, char* argv[]) {

  high_resolution_clock::time_point start;

  // used bool to save space, denoting is_director
  vector<tuple<unsigned int, unsigned int, unsigned short> > principals = load_principals("/bigstore/yike/yqiuac/imdb/principals.tsv");
  vector<tuple<unsigned int, double> > ratings = load_ratings("/bigstore/yike/yqiuac/imdb/ratings.tsv");

  cerr << "-- principals ---" << endl;
  unsigned int N_pri = principals.size();
  cerr << "INPUT SIZE: " << N_pri << endl;
  unsigned int D_pri = get_distinct_count(principals);
  cerr << "DISTINCT COUNT: " << D_pri << endl;

  
  // ////////////////////////////// EXACT JOIN //////////////////////////////

  set_timer(start);
  
  RELATION_TYPE data; // abcxyz

  sort(ratings.begin(), ratings.end());
  for (const auto& principals_t : principals) {
    auto nconst = get<0>(principals_t);
    auto tconst = get<1>(principals_t);
    auto category = get<2>(principals_t);
    auto it_begin = lower_bound
      (ratings.begin(), ratings.end(), tconst,
       &tuple_less_key<unsigned int, double>);  // quivalent to find in sorted list
    if (it_begin == ratings.end()) {
      continue;
    }
    auto it_end = upper_bound
      (ratings.begin(), ratings.end(), tconst,
       &key_less_tuple<unsigned int, double>);  // quivalent to find in sorted list

    for (auto it = it_begin; it != it_end; ++it) {
      auto rating = get<1>(*it);
      data.emplace_back(make_tuple(nconst, tconst, category, rating));
    }
  }

  auto time_span_join = get_timer(start);

  cerr << "-- joint table ---" << endl;
  unsigned int N = data.size();
  cerr << "INPUT SIZE: " << N << endl;
  unsigned int D = get_distinct_count(data);
  cerr << "DISTINCT COUNT: " << D << endl;
  
  cerr << "-- sample ---" << endl;
  unsigned int n = round(0.01 * N);
  cerr << "SAMPLE SIZE: " << n << endl;

  // Solve for WDS parameters
  unsigned int M; // Number of 1's
  unordered_map<unsigned int, double> p_i; // Prob of each value
  duration<double> time_span_solve = solve(data, D, n, M, p_i);

  unsigned int M_rw; // Number of 1's
  unordered_map<unsigned int, double> p_i_rw; // Prob of each value
  duration<double> time_span_solve_rw = solve(data, D, n, M_rw, p_i_rw);

  /**
   * Each file will contain 30 lines, each line is an independent trial
   * Each trial will contain 9 values for the 9 filters
   * year moy class category (*30)
   */
  ofstream uds_file;
  uds_file.open("../result/uds.txt");
  ofstream wds_file;
  wds_file.open("../result/wds.txt");
  ofstream exact_file;
  exact_file.open("../result/exact.txt");
  ofstream ub_file;
  ub_file.open("../result/ub.txt");
  ofstream rw_file;
  rw_file.open("../result/rw.txt");

  /**
   * time_file will contain (1+9)*30 lines
   * Every 10 lines, the first line is the sample construction time for uds/wds/rw
   * The remaining 9 lines is the query time for uds/wds/exact/rw
   */

  ofstream time_file; 
  time_file.open("../result/time.txt");
  
  for (int trial = 1; trial <= 30; ++trial) {
    cerr << "Running Trial " << trial;
  
    srand(time(NULL));
  
    ////////////////////////////// UDS //////////////////////////////
  
    set_timer(start);
  
    unsigned int t_ds = round(n / 50.0);
    unsigned int l_ds; 
    RELATION_TYPE uds_sketch = collect_uds_sketch(data, n, t_ds, l_ds);

    cerr << '.';
    auto time_span_uds = time_span_join + get_timer(start);

    ////////////////////////////// WDS //////////////////////////////

    set_timer(start);
    
    unordered_map<unsigned int, RELATION_TYPE> wds_sketch = collect_wds_sketch(data, p_i);

    cerr << '.';
    auto time_span_wds = time_span_join + time_span_solve + get_timer(start);

    ////////////////////////////// Random Walk //////////////////////////////

    set_timer(start);

    // Process sales table
    // get a sample of items
    // and record their frequency count (which will be the upperbound)

    sort(ratings.begin(), ratings.end());

    double_hasher h;
    unordered_map<unsigned int, sketch_unit> rw_sketch;

    // First Table
    for (const auto& principals_t : principals) {
      auto nconst = get<0>(principals_t);

      if (h(nconst) <= p_i_rw[nconst]) { // This value should be sampled
  	rw_sketch[nconst].add_table(principals_t);
      }
    }

    for (auto& sketch_pair : rw_sketch) {
      sketch_pair.second.commit_table();
    }

    // Second table
    for (auto& sketch_pair : rw_sketch) { 
      auto& sketch_a = sketch_pair.second;

      for (const auto& table_t : sketch_a.table) {
  	auto tconst = get<1>(table_t);

  	auto it_begin = lower_bound
  	  (ratings.begin(), ratings.end(), tconst,
  	   &tuple_less_key<unsigned int, double>);  // quivalent to find in sorted list
  	if (it_begin == ratings.end()) {
  	  continue;
  	}
  	auto it_end = upper_bound
  	  (ratings.begin(), ratings.end(), tconst,
  	   &key_less_tuple<unsigned int, double>);  // quivalent to find in sorted list

  	for (auto it = it_begin; it != it_end; ++ it) {
  	  sketch_a.add_binary_join(table_t, *it);
  	}

      }
      sketch_a.commit_binary_join();
    }

    cerr << '.';
    auto time_span_rw = time_span_solve_rw + get_timer(start);

    time_file << time_span_uds.count() << '\t'
    	      << time_span_wds.count() << '\t'
    	      << time_span_rw.count() << endl;
    cerr << endl;

    ////////////////////////////// ESTIMATE //////////////////////////////

    for (int param = 1; param <= 9; ++param) {

      ////////////////////////////// UDS Estimate //////////////////////////////

      set_timer(start);

      RELATION_TYPE uds_sketch_filtered = filter_table(uds_sketch, param);
      unsigned int uds_sketch_distinct_count = get_distinct_count(uds_sketch_filtered);
      double uds_estimator = uds_sketch_distinct_count * pow(2, l_ds);
      uds_sketch_filtered.clear();

      time_span_uds = get_timer(start);
  
      ////////////////////////////// WDS Estimate //////////////////////////////

      set_timer(start);
      
      double wds_estimator = 0.0;
      for (auto it = wds_sketch.cbegin(); it != wds_sketch.cend(); ++it) {
  	RELATION_TYPE Ra = it -> second;
  	RELATION_TYPE Ra_filtered = filter_table(Ra, param);
  	if (Ra_filtered.size() > 0) {
  	  wds_estimator += 1.0 / p_i[it -> first];
  	}
  	Ra_filtered.clear();
      }
      //      wds_estimator = wds_estimator + (D - M) / 2.0;

      time_span_wds = get_timer(start);

      ////////////////////////////// Exact //////////////////////////////

      set_timer(start);

      RELATION_TYPE filtered_data = filter_table(data, param);
      unsigned int filtered_size = filtered_data.size();
      unsigned int filtered_distinct_count = get_distinct_count(filtered_data);
      filtered_data.clear();

      duration<double> time_span_exact = time_span_join + get_timer(start);

      ////////////////////////////// Random Walk Estimate //////////////////////////////

      set_timer(start);

      double rw_estimator = 0.0;
      for (const auto& sketch_pair : rw_sketch) {
  	const auto& a = sketch_pair.first;
  	const auto& sketch_a = sketch_pair.second;

  	rw_estimator += sketch_a.estimate(param) / p_i_rw[a];
      }
      //rw_estimator = rw_estimator + (D_pri - M_rw) / 2.0;
      time_span_rw = get_timer(start);

      ////////////////////////////// UB //////////////////////////////

      double ub_estimator = min(filtered_size, D_pri);

      ////////////////////////////// OUTPUT //////////////////////////////

      time_file << time_span_uds.count() << '\t'
  		<< time_span_wds.count() << '\t'
  		<< time_span_exact.count() << '\t'
  		<< time_span_rw.count() << endl;

      exact_file << filtered_distinct_count << '\t';
      ub_file << ub_estimator << '\t';
      uds_file << uds_estimator << '\t';
      wds_file << wds_estimator << '\t';
      rw_file << rw_estimator << '\t';
    }
    exact_file << endl;
    ub_file << endl;
    uds_file << endl;
    wds_file << endl;
    rw_file << endl;
    cerr << endl;
  }
  exact_file.close();
  ub_file.close();
  uds_file.close();
  wds_file.close();
  rw_file.close();
  time_file.close();

  return 0;
}
