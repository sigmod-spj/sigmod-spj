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

typedef tuple<unsigned int, unsigned int>  TUPLE_TYPE;
typedef tuple<unsigned int, unsigned int, unsigned int>  TUPLE_12_TYPE;
typedef tuple<unsigned int, unsigned int, unsigned int, unsigned int>  TUPLE_123_TYPE;
typedef tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>  TUPLE_1234_TYPE;
typedef tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>  TUPLE_12345_TYPE;

typedef vector<tuple<unsigned int, unsigned int> > TABLE_TYPE;
typedef vector<tuple<unsigned int, unsigned int, unsigned int> > TABLE_12_TYPE;
typedef vector<tuple<unsigned int, unsigned int, unsigned int, unsigned int> > TABLE_123_TYPE;
typedef vector<tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int> > TABLE_1234_TYPE;
typedef vector<tuple<unsigned int, unsigned int, unsigned int, unsigned int > > RELATION_TYPE;

/* Q54 Filter From Agirish */
template <class K, class... V>
bool tuple_filter(const tuple<K, V...>& t, int param) {
  auto& a = get<0>(t);
  auto& b = get<1>(t);
  auto& c = get<2>(t);
  auto& d = get<3>(t);
  return (d==1) && (b != a) &&
    (c != a) && (c != b) &&
    (d != a) && (d != b) && (d != c);
}

struct sketch_unit {
  unsigned int max_freq;
  unsigned int current_freq;
  bool was_complete;

  vector<double> prob;
  vector<double> temp_prob;

  TABLE_TYPE table_1;
  TABLE_12_TYPE table_12;
  TABLE_123_TYPE table_123;

  sketch_unit(): max_freq(0), current_freq(0), was_complete(false) {}
  ~sketch_unit() {
    prob.clear();
    temp_prob.clear();
    table_1.clear();
    table_12.clear();
    table_123.clear();
  }

  double estimate(int param) const {
    bool success = false;
    double p = 1.0;
    for (unsigned int i = 0; i < table_123.size(); ++i) {
      if (tuple_filter(table_123[i], param)) {
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

  unsigned int commit_table_1() {
    current_freq = 0;
    //    max_freq = max_freq * 2;
    return table_1.size();
  }
  unsigned int commit_table_12() {
    if (!was_complete && current_freq > max_freq) {
      was_complete = true;
      prob = vector<double>(table_12.size(), max_freq / static_cast<double>(current_freq));
    }
    table_1.clear();
    current_freq = 0;
    return table_12.size();
  }
  unsigned int commit_table_123() {
    if (was_complete) {
      prob = temp_prob;
      temp_prob.clear();
    } else if (current_freq > max_freq) {
      was_complete = true;
      prob = vector<double>(table_123.size(), max_freq / static_cast<double>(current_freq));
    }
    table_12.clear();
    current_freq = 0;
    return table_123.size();
  }


  void add_table_1(const TUPLE_TYPE &s) {
    table_1.emplace_back(s);
    current_freq += 1;
    max_freq += 1;
  }

  // abcxy
  void add_table_12(const TUPLE_TYPE &s,
                    const TUPLE_TYPE &d) {
    if (current_freq < max_freq) {
      table_12.emplace_back
        (make_tuple(get<0>(s), get<1>(s), get<1>(d)));
    } else {
      unsigned int idx = rand() % current_freq;
      if (idx < max_freq) {
    	table_12[idx] = make_tuple(get<0>(s), get<1>(s), get<1>(d));
      }
    }
    current_freq += 1;
  }

  void add_table_123(const TUPLE_12_TYPE &sd,
                     const TUPLE_TYPE &i) {
    if (current_freq < max_freq) {
      table_123.emplace_back
        (make_tuple(get<0>(sd), get<1>(sd), get<2>(sd), get<1>(i)));
    } else {
      unsigned int idx = rand() % current_freq;
      if (idx < max_freq) {
        table_123[idx] = make_tuple(get<0>(sd), get<1>(sd), get<2>(sd), get<1>(i));
      }
    }
    current_freq += 1;
  }
};

ostream& operator<<(ostream &os, sketch_unit s) {
  os << "MAX FREQUENCY: " << s.max_freq << endl
     << "CURRENT FREQUENCY: " << s.current_freq << endl
     << "WAS COMPLETE: " << s.was_complete << endl
     << "TABLE SIZE: " << s.table_1.size() << endl
     << "TABLE 12 SIZE: " << s.table_12.size() << endl
     << "TABLE 123 SIZE: " << s.table_123.size() << endl
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

// solve will sort the data
template <class KK, class... V>
void solve
(vector<tuple<KK, V...> > &data,
 double D,
 double n,
 unsigned int &M,
 unordered_map<KK, double> &p_i) {
  if (data.empty()) {
    return;
  }
  sort(data.begin(), data.end());

  // Collect Degree Information
  vector<pair<unsigned int, KK> > N_i;
  N_i.emplace_back(make_pair(0, KK())); // template for N_i[0], so that we can directly call N_i[i].

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

  double min_error = std::numeric_limits<double>::max();
  double min_M = 0;
  double min_K = 0;

  double error = D * D;
  unsigned int M_0 = 0;

  double sum = 0.0;

  for (auto i = 1; i <= D; ++i) {
    sum += N_i[i].first;
    if (n >= sum) {
      M_0 = i;
    } else {
      sum -= N_i[i].first;
      break;
    }
  } // EndFor
  // Now M_0 is the largest integer that n >= sum_1^M0 N_i

  unsigned int K = M_0; // sum is sum_1^K N_i
  double N_k = N_i[K].first;
  double sum_sqrt = -sqrt(N_k); // sum_sqrt will be sum_K+1^M sqrt(N_i)


  for (M = M_0; M <= D; ++M) {
    sum_sqrt += sqrt(N_i[M].first);
    while (K > 0 && n - sum <= sqrt(N_k) * sum_sqrt) { // K is too large
      sum -= N_k;
      sum_sqrt += sqrt(N_k);
      K = K-1;
      N_k = N_i[K].first;
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

  for (int i = 1; i <= K; ++i) {
    sum += N_i[i].first;
    p_i.emplace(N_i[i].second, 1.0);
  } // 1...K

  for (int i = K+1; i <= M; ++i) {
    sum_sqrt += sqrt(N_i[i].first);
  }
  for (int i = K+1; i <= M; ++i) {
    p_i.emplace(N_i[i].second, (n - sum) / sum_sqrt / sqrt(N_i[i].first));
  } // K+1 ... M

  for (int i = M+1; i <= D; ++i) {
    p_i.emplace(N_i[i].second, 0.0);
  } // M+1 ... D

  cerr << "Parameter K = " << K << ", M = " << M << endl;
}


template <class K, class... V>
vector<tuple<K, V...> > collect_ds_sketch
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
unordered_map<K, vector<tuple<K, V...> > > collect_dads_sketch
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

  using namespace std::chrono;
  high_resolution_clock::time_point start;

  TABLE_TYPE graph = load_graph((argc>1)?argv[1]:"../data/real/p2p-Gnutella04.txt");

  cerr << "-- graph ---" << endl;
  unsigned int N_graph = graph.size();
  cerr << "INPUT SIZE: " << N_graph << endl;
  unsigned int D_graph = get_distinct_count(graph);
  cerr << "DISTINCT COUNT: " << D_graph << endl;


  // ////////////////////////////// EXACT JOIN //////////////////////////////

  set_timer(start);

  RELATION_TYPE data; // abcxyz

  TABLE_12_TYPE table_12;
  TABLE_123_TYPE table_123;

  for (const auto& graph_ab : graph) {
    auto a = get<0>(graph_ab);
    auto b = get<1>(graph_ab);
    auto b_it_begin = lower_bound
      (graph.begin(), graph.end(), b,
       &tuple_less_key<unsigned int, unsigned int>);  // quivalent to find in sorted list
    if (b_it_begin == graph.end()) {
      continue;
    }
    auto b_it_end = upper_bound
      (graph.begin(), graph.end(), b,
       &key_less_tuple<unsigned int, unsigned int>);  // quivalent to find in sorted list

    for (auto b_it = b_it_begin; b_it != b_it_end; ++b_it) {
      auto c = get<1>(*b_it);
      table_12.emplace_back(make_tuple(a, b, c));
    }
  }

  for (const auto& graph_abc : table_12) {
    auto a = get<0>(graph_abc);
    auto b = get<1>(graph_abc);
    auto c = get<2>(graph_abc);
    auto c_it_begin = lower_bound
      (graph.begin(), graph.end(), c,
       &tuple_less_key<unsigned int, unsigned int>);  // quivalent to find in sorted list
    if (c_it_begin == graph.end()) {
      continue;
    }
    auto c_it_end = upper_bound
      (graph.begin(), graph.end(), c,
       &key_less_tuple<unsigned int, unsigned int>);  // quivalent to find in sorted list

    for (auto c_it = c_it_begin; c_it != c_it_end; ++c_it) {
      auto d = get<1>(*c_it);
      data.emplace_back(make_tuple(a, b, c, d));
    }
  }
  table_12.clear();


  auto time_span_join = get_timer(start);

  cerr << "-- joint table ---" << endl;
  unsigned int N = data.size();
  cerr << "INPUT SIZE: " << N << endl;
  unsigned int D = get_distinct_count(data);
  cerr << "DISTINCT COUNT: " << D << endl;

  cerr << "-- sketch ---" << endl;
  unsigned int n = round(0.01 * N);
  cerr << "SAMPLE SIZE: " << n << endl;

  unsigned int M; // Number of 1's
  unordered_map<unsigned int, double> p_i; // Prob of each value
  solve(data, D, n, M, p_i);

  unsigned int M_rw; // Number of 1's
  unordered_map<unsigned int, double> p_i_rw; // Prob of each value
  solve(graph, D_graph, n, M_rw, p_i_rw);



  cerr << "EXACT" << '\t'
       << "UB" << '\t'
       << "DS" << '\t'
       << "DADS" << '\t'
       << "RW" <<'\t'<<endl;


  /**
   * Format: (*trial)
   * filter 1 / 2 / ... / param
   */
  ofstream ds_file;
  ds_file.open("../result/ds.txt");
  ofstream dads_file;
  dads_file.open("../result/dads.txt");
  ofstream exact_file;
  exact_file.open("../result/exact.txt");
  ofstream ub_file;
  ub_file.open("../result/ub.txt");
  ofstream rw_file;
  rw_file.open("../result/rw.txt");


  /**
   * Format: (*trial)
   * Unit: s
   * DS_SKETCH_TIME / RP_SKETCH_TIME / DADS_SKETCH_TIME / JOIN_TIME / RW_SKETCH_TIME (1)
   * DS_QUERY_TIME / RP_QUERY_TIME / DADS_QUERY_TIME / EXACT_QUERY_TIME (*param)
   */
  ofstream time_file;
  time_file.open("../result/time.txt");

  for (int trial = 1; trial <= 30; ++trial) {
    cerr << "Running Trial " << trial;

     srand(time(NULL));

    ////////////////////////////// DS //////////////////////////////

    set_timer(start);

    unsigned int t_ds = round(n / 50.0);
    unsigned int l_ds;
    RELATION_TYPE ds_sketch = collect_ds_sketch(data, n, t_ds, l_ds);

    cerr << '.';
    auto time_span_ds = get_timer(start);

    ////////////////////////////// DADS //////////////////////////////

    set_timer(start);

    unordered_map<unsigned int, RELATION_TYPE> dads_sketch = collect_dads_sketch(data, p_i);

    cerr << '.';
    auto time_span_dads = get_timer(start);

    ////////////////////////////// Random Walk //////////////////////////////

    set_timer(start);

    // Process sales table
    // get a sample of items
    // and record their frequency count (which will be the upperbound)

    sort(graph.begin(), graph.end());

    double_hasher h;
    unordered_map<unsigned int, sketch_unit> rw_sketch;

    unsigned int rw_size = 0;
    // First Table
    for (const auto& graph_t : graph) {
      auto a = get<0>(graph_t);

      if (h(a) <= p_i_rw[a]) { // This value should be sampled
        rw_sketch[a].add_table_1(graph_t);
      }
    }

    for (auto& sketch_pair : rw_sketch) {
      rw_size += sketch_pair.second.commit_table_1();
    }

    //    cerr << "RW sketch size of table 1: " << rw_size << " / n: " << n << endl;
    rw_size = 0;


    // Second table
    for (auto& sketch_pair : rw_sketch) {
      auto& sketch_a = sketch_pair.second;

      for (const auto& table_t : sketch_a.table_1) {
        auto b = get<1>(table_t);

        auto b_it_begin = lower_bound
          (graph.begin(), graph.end(), b,
           &tuple_less_key<unsigned int, unsigned int>);
        if (b_it_begin == graph.end()) {
          continue;
        }
        auto b_it_end = upper_bound
          (graph.begin(), graph.end(), b,
           &key_less_tuple<unsigned int, unsigned int>);

        for (auto b_it = b_it_begin; b_it != b_it_end; ++ b_it) {
          sketch_a.add_table_12(table_t, *b_it);
        }
      }
      rw_size += sketch_a.commit_table_12();
    }

    //    cerr << "RW sketch size of table 12: " << rw_size << " / n: " << n << endl;
    rw_size = 0;

    //3rd
    for (auto& sketch_pair : rw_sketch) {
      auto& sketch_a = sketch_pair.second;

      for (unsigned int i = 0; i < sketch_a.table_12.size(); ++i) {
        const auto& sketch_t = sketch_a.table_12[i];
        auto c = get<2>(sketch_t);

        auto c_it_begin = lower_bound
          (graph.begin(), graph.end(), c,
           &tuple_less_key<unsigned int, unsigned int>);
        if (c_it_begin == graph.end()) {
          continue;
        }
        auto c_it_end = upper_bound
          (graph.begin(), graph.end(), c,
           &key_less_tuple<unsigned int, unsigned int>);
        if (c_it_begin == c_it_end) {
          continue;
        }
        if (sketch_a.was_complete) { // RANDOM WALK STAGE
          auto tuple_prob = sketch_a.prob[i];
          auto fan_out = c_it_end - c_it_begin;
          unsigned int idx = rand() % fan_out;
          auto c_it = c_it_begin + idx;
          sketch_a.add_table_123(sketch_t, *c_it);
          sketch_a.temp_prob.emplace_back(tuple_prob / fan_out);
        } else {
          for (auto c_it = c_it_begin; c_it != c_it_end; ++ c_it) {
            sketch_a.add_table_123(sketch_t, *c_it);
          } // END FOR ITEM_IT
        } // ENDIF
      } // END FOR SALES
      rw_size += sketch_a.commit_table_123();
    } // END FOR SKETCH_IT



    cerr << "RW sketch size of table 12345: " << rw_size << " / n: " << n << endl;

    // int print_count = 0;
    // for (auto iter = rw_sketch.begin(); iter!=rw_sketch.end() && print_count < 10; ++iter, ++print_count) {
    //   cout << iter -> first << endl << iter->second << endl;
    // }


    auto time_span_rw = get_timer(start);

    time_file << time_span_ds.count() << '\t'
              << time_span_dads.count() << '\t'
              << time_span_join.count() << '\t'
              << time_span_rw.count()
              << endl;


    ////////////////////////////// ESTIMATE //////////////////////////////

    int param = 1;

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();

    ////////////////////////////// DS Estimate //////////////////////////////

    set_timer(start);

    RELATION_TYPE ds_sketch_filtered = filter_table(ds_sketch, param);
    unsigned int ds_sketch_distinct_count = get_distinct_count(ds_sketch_filtered);
    double ds_estimator = ds_sketch_distinct_count * pow(2, l_ds);
    ds_sketch_filtered.clear();

    time_span_ds = get_timer(start);


      ////////////////////////////// DADS Estimate //////////////////////////////

    set_timer(start);

    double dads_estimator = 0.0;
    for (auto it = dads_sketch.cbegin(); it != dads_sketch.cend(); ++it) {
      RELATION_TYPE Ra = it -> second;
      RELATION_TYPE Ra_filtered = filter_table(Ra, param);
      if (Ra_filtered.size() > 0) {
        dads_estimator += 1.0 / p_i[it -> first];
      }
      Ra_filtered.clear();
    }
    // dads_estimator += (D*1.0-M) / 2.0;

    time_span_dads = get_timer(start);

    ////////////////////////////// Exact //////////////////////////////

    set_timer(start);

    RELATION_TYPE filtered_data = filter_table(data, param);
    unsigned int filtered_size = filtered_data.size();
    unsigned int filtered_distinct_count = get_distinct_count(filtered_data);
    filtered_data.clear();

    duration<double> time_span_exact = get_timer(start);

    ////////////////////////////// Random Walk Estimate //////////////////////////////

    set_timer(start);

    double rw_estimator = 0.0;
    for (const auto& sketch_pair : rw_sketch) {
      const auto& a = sketch_pair.first;
      const auto& sketch_a = sketch_pair.second;

      rw_estimator += sketch_a.estimate(param) / p_i_rw[a];
    }
    rw_estimator = rw_estimator + (D_graph*1.0 - M_rw) / 2.0;

    time_span_rw = get_timer(start);

    ////////////////////////////// UB //////////////////////////////

    double ub_estimator = min(filtered_size, D_graph);

    ////////////////////////////// OUTPUT //////////////////////////////

    time_file << time_span_ds.count() << '\t'
              << time_span_dads.count() << '\t'
              << time_span_exact.count() << '\t'
              << time_span_rw.count() << endl;

    cerr << filtered_distinct_count << '\t'
         << ub_estimator << '\t'
         << ds_estimator << '\t'
         << dads_estimator << '\t'
         << rw_estimator;


    exact_file << filtered_distinct_count << '\t';
    ub_file << ub_estimator << '\t';
    ds_file << ds_estimator << '\t';
    dads_file << dads_estimator << '\t';
    rw_file << rw_estimator << '\t';

    exact_file << endl;
    ub_file << endl;
    ds_file << endl;
    dads_file << endl;
    rw_file << endl;
    cerr << endl;
   }
  exact_file.close();
  ub_file.close();
  ds_file.close();
  dads_file.close();
  rw_file.close();
  time_file.close();

  return 0;
}
