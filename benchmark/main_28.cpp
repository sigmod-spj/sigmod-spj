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

// There are 6 different filters in each trial
#define PARAM_UB 6

// Type of attributes we take from the store_sales table
typedef std::vector<std::tuple<unsigned int, unsigned int, double, double> > RELATION_TYPE;

using namespace std;
using namespace std::chrono;

/**
 * TPC-DS Q28 Filter
 * @param t input tuple
 * @param i Bucket number B1-B6 for i=0-5
 * @return whether tuple t passes the filter
 */

template <class K, class... V>
bool tuple_filter(const tuple<K, V...>& t, unsigned int i) {
  double ss_list_price = get<0>(t) / 100.0;
  unsigned int ss_quantity = get<1>(t);
  double ss_coupon_amt = get<2>(t);
  double ss_wholesale_cost = get<3>(t);

  unsigned int ss_quantity_low;
  unsigned int ss_list_price_low;
  unsigned int ss_coupon_amt_low;
  unsigned int ss_wholesale_cost_low;

  switch (i) {
  case 0:
    ss_quantity_low = 1;
    ss_list_price_low = 18;
    ss_coupon_amt_low = 1939;
    ss_wholesale_cost_low = 34;
    break;
  case 1:
    ss_quantity_low = 6;
    ss_list_price_low = 1;
    ss_coupon_amt_low = 35;
    ss_wholesale_cost_low = 50;
    break;
  case 2:
    ss_quantity_low = 11;
    ss_list_price_low = 91;
    ss_coupon_amt_low = 1412;
    ss_wholesale_cost_low = 17;
    break;
  case 3:
    ss_quantity_low = 16;
    ss_list_price_low = 9;
    ss_coupon_amt_low = 5270;
    ss_wholesale_cost_low = 29;
    break;
  case 4:
    ss_quantity_low = 21;
    ss_list_price_low = 45;
    ss_coupon_amt_low = 826;
    ss_wholesale_cost_low = 5;
    break;
  case 5:
    ss_quantity_low = 26;
    ss_list_price_low = 174;
    ss_coupon_amt_low = 5548;
    ss_wholesale_cost_low = 42;
    break;
  }
  return (ss_quantity >= ss_quantity_low) &&
    (ss_quantity <= ss_quantity_low + 4) &&
    ((ss_list_price >= ss_list_price_low && ss_list_price <= ss_list_price_low + 10) ||
     (ss_coupon_amt >= ss_coupon_amt_low && ss_coupon_amt <= ss_coupon_amt_low + 1000) ||
     (ss_wholesale_cost >= ss_wholesale_cost_low && ss_wholesale_cost <= ss_wholesale_cost_low + 20) );
}

/**
 * generates a filtered table
 * @param table the input table
 * @param i parameter for the filter
 * @return a filtered table that only contains tuples in table that passes the filter of i.
 */
RELATION_TYPE filter_table(const RELATION_TYPE& table,unsigned int i) {
  RELATION_TYPE result;
  for (auto it = table.cbegin(); it != table.cend(); ++it) {
    if (tuple_filter(*it,i)) {
      result.emplace_back(*it);
    }
  }
  return result;
}

/**
 * Sort the table and then count the distinct keys in one pass
 * @param data the related to be counted, will be sorted in-place
 * @return the NDV of key of data
 */
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

/**
 * The histogram algorithm described in section 4.2
 */
template <class KK, class... V>
duration<double> solve_histo
(vector<tuple<KK, V...> > &data,
 double D,
 unsigned int n,
 unsigned int &M,
 unordered_map<KK, double> &p_i,
 double &p_low) {
  if (data.empty()) {
    return duration<double>();
  }
  sort(data.begin(), data.end());

  double sum_freq = 0.0;
  // Collect Degree Information
  vector<pair<double, KK> > N_i;
  unsigned int N_a = 1;
  KK current(get<0>(*data.cbegin()));
  for (auto it = ++data.cbegin(); it != data.cend(); ++it) {
    if (current == get<0>(*it)) {
      N_a += 1;
    } else {
      sum_freq += N_a;
      N_i.emplace_back(make_pair(N_a, current));
      current = get<0>(*it);
      N_a = 1;
    }
  }
  sum_freq += N_a;
  N_i.emplace_back(make_pair(N_a, current));

  // Sort value by degrees
  sort(N_i.begin(), N_i.end()); // increasing


  double sum_heavy = 0.0;
  unsigned int H = static_cast<int>(D * 9.0 /10.0);
  for (unsigned int idx = H; idx < D-0.5; ++idx) { // N_H + ... + N_{D-1}
    sum_heavy += N_i[idx].first;
  }
  double sum_light = sum_freq - sum_heavy; // N_0 + ... + N_H-1
  double N_avg = sum_light / H; // N_0 ... N_{H-1}

  high_resolution_clock::time_point start = high_resolution_clock::now();

  double min_error = std::numeric_limits<double>::max();
  double min_M = 0;
  double min_K = 0;

  double error = D * D;
  unsigned int M_0 = 0;

  double sum = 0.0;
  unsigned int K;

  if (sum_light <= n) { // N_1 + ... + N_H <= n, so M_0 >= H
    sum = sum_light;
    for (M_0 = H+1; M_0 <= D; ++M_0) {
      sum += N_i[M_0-1].first; // sum += N_H ... sum += N_{D-1}
      if (sum > n) {
        sum -= N_i[M_0-1].first;
        M_0 = M_0 - 1;
        break;
      }
    }
  } else { // M_0 is in 1 ... H
    M_0 = floor(n / N_avg);
    sum = M_0 * N_avg;

    cout << "M_0 = " << M_0 << endl;

    // Solution 1: M<=H and n <= M * N_avg
    M = round((2*D + 1) / (2 * (1 + N_avg / n)));
    cout << "Solution 1: M = " << M <<endl;

    if (M<=H && n <= M*N_avg) {
      K = 0;
      error = (D-M) * (D-M) + N_avg / n * M * M  -M;
      if (error < min_error) {
        min_error = error;
        min_M = M;
        min_K = K;
      }
    }
    // Solution 2: M<=H and n > M * N_avg
    M = floor(n / N_avg);
    cout << "Solution 2: M = " << M <<endl;
    if (M<=H && n > M*N_avg) {
      K = M;
      error = (D-M) * (D-M);
      if (error < min_error) {
        min_error = error;
        min_M = M;
        min_K = K;
      }
    }

  }
  // Now M_0 is the largest integer that n >= sum_1^M0 N_i
  // sum is N_1 + ... + N_{M_0}
  // We only need to consider M_0 > H

  // Case 2: M > H
  // Only need to search M >= H and K <= M_0
  K = M_0; // sum is sum_1^K N_i
  double N_k = N_i[K-1].first;
  double sum_sqrt = 0.0;
  if (M_0 >= H) {
    M = M_0;
    sum_sqrt = -sqrt(N_k); // sum_sqrt will be sum_K+1^M sqrt(N_i)
  } else {
    M = H+1;
    sum_sqrt = (sum_light - sum) / sqrt(N_avg);
  }

  cout << "sum: " << sum << endl;
  cout << "sum_sqrt: " << sum_sqrt << endl;

  for (; M <= D; ++M) {
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

  if (K >= H) {
    p_low = 1.0;
  } else if (K == 0) {
    p_low = 0.0;
  }
  for (unsigned int i = H; i < K; ++i) {
    sum += N_i[i].first;
    p_i.emplace(N_i[i].second, 1.0);
  } // 1...K

  for (unsigned int i = K; i < M; ++i) {
    sum_sqrt += sqrt(N_i[i].first);
  }
  for (unsigned int i = K; i < M; ++i) {
    p_i.emplace(N_i[i].second, (n - sum) / sum_sqrt / sqrt(N_i[i].first));
  } // K+1 ... M

  for (unsigned int i = M; i < D; ++i) {
    p_i.emplace(N_i[i].second, 0.0);
  }

  cerr << "Parameter K_histo = " << K << ", M_histo = " << M << endl;
  high_resolution_clock::time_point end = high_resolution_clock::now();
  return duration_cast<duration<double>>(end - start);

}

/**
 * Collect a uniform distinct sampling sketch
 * @return the sketch as a relation
 */
RELATION_TYPE collect_uds_sketch
(const RELATION_TYPE &data,
 unsigned int n, unsigned int t, unsigned int &l) {
  unordered_map<unsigned int, RELATION_TYPE> S;
  unsigned int S_size = 0;
  unordered_map<unsigned int, unsigned int> c;
  die_hasher h;
  l = 0;
  for (auto it = data.cbegin(); it != data.cend(); ++ it) {
    unsigned int v = get<0>(*it);
    unsigned int die_level = h(v);
    if (die_level >= l) {
      if (S.find(v)==S.end()) {
        S.emplace(v, RELATION_TYPE());
        S[v].emplace_back(*it);
        c[v] = 1;
        S_size += 1;
      } else if (c[v] <= t-1) {
        S[v].emplace_back(*it);
        c[v] += 1;
        S_size += 1;
      } else {
        unsigned int idx = rand() % c[v];
        if (idx < t) {
          S[v][idx] = *it;
        }
        c[v] += 1;
      }

      if (S_size >= n) {
        for (auto map_it = S.begin(); map_it != S.end();) {
          unsigned int v = map_it -> first;
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
  RELATION_TYPE result;
  for (auto map_it = S.begin(); map_it != S.end(); ++map_it) {
    result.insert(result.end(),map_it->second.begin(),map_it->second.end());
  }
  return result;
}

/**
 * Collect the weighted distinct sampling sketch similar to example 2
 */
unordered_map<unsigned int, RELATION_TYPE> collect_wds_sketch
(const RELATION_TYPE &data,
 const unordered_map<unsigned int,double>& p_i) {
  unordered_map<unsigned int, RELATION_TYPE> S;
  unsigned int S_size = 0;
  unordered_map<unsigned int, unsigned int> c;
  double_hasher h;

  for (auto it = data.cbegin(); it != data.cend(); ++ it) {
    unsigned int a = get<0>(*it);
    double v = h(a);
    if (v <= p_i.find(a) -> second) {
      if (S.find(a)!=S.end()) {
        S[a].emplace_back(*it);
        c[a] += 1;
        S_size += 1;
      } else {
        S.emplace(a, RELATION_TYPE());
        S[a].emplace_back(*it);
        c[a] = 1;
        S_size += 1;
      }
    }
  }
  return S;
}

int main(int argc, char* argv[]) {

  RELATION_TYPE data = load_store_sales((argc>1)?argv[1]:"../data/tpcds/store_sales.dat");

  cout << "-- data ---" << endl;
  unsigned int N = data.size();
  cout << "INPUT SIZE: " << N << endl;
  unsigned int D = get_distinct_count(data);
  cout << "DISTINCT COUNT: " << D << endl;

  cout << "-- sample ---" << endl;
  unsigned int n = round(0.01 * N);
  cout << "SAMPLE SIZE: " << n << endl;

  // Solve for WDS parameters
  unsigned int M; // Number of 1's
  unordered_map<unsigned int, double> p_i; // Prob of each value
  duration<double> time_span_solve = solve(data, D, n, M, p_i);

  unsigned int M_histo; // Number of 1's
  unordered_map<unsigned int, double> p_i_histo; // Prob of each value
  double p_low;
  duration<double> time_span_solve_histo = solve_histo(data, D, n, M_histo, p_i_histo, p_low);

  /**
   * Each file will contain 30 lines, each line is an independent trial
   * Each trial will contain 6 values for the 6 filters
   * B1 B2 B3 B4 B5 B6 (*30)
   */
  ofstream uds_file;
  uds_file.open("../result/uds.txt");
  ofstream wds_file;
  wds_file.open("../result/wds.txt");
  ofstream exact_file;
  exact_file.open("../result/exact.txt");
  ofstream ub_file;
  ub_file.open("../result/ub.txt");
  ofstream histo_file;
  histo_file.open("../result/histo.txt");

  /**
   * time_file will contain (1+6)*30 lines
   * Every 7 lines, the first line is the sample construction time for uds/wds/histo
   * The remaining 6 lines is the query time for uds/wds/histo/exact
   */
  ofstream time_file;
  time_file.open("../result/time.txt");

  for (int trial = 1; trial <= 30; ++trial) {
    cout << "Running Trial " << trial;

    srand(time(NULL));

    ////////////////////////////// UDS //////////////////////////////

    high_resolution_clock::time_point start = high_resolution_clock::now();

    unsigned int t_ds = round(n / 50.0);
    unsigned int l_ds;
    RELATION_TYPE uds_sketch = collect_uds_sketch(data, n, t_ds, l_ds);

    high_resolution_clock::time_point end = high_resolution_clock::now();
    cout << '.';
    duration<double> time_span_uds = duration_cast<duration<double>>(end - start);

    ////////////////////////////// WDS //////////////////////////////

    start = high_resolution_clock::now();

    unordered_map<unsigned int, RELATION_TYPE> wds_sketch = collect_wds_sketch(data, p_i);

    end = high_resolution_clock::now();
    cout << '.';
    duration<double> time_span_wds = time_span_solve + duration_cast<duration<double>>(end - start);

    ////////////////////////////// HISTO //////////////////////////////

    start = high_resolution_clock::now();

    unordered_map<unsigned int, RELATION_TYPE> histo_sketch = collect_wds_sketch(data, p_i_histo);

    end = high_resolution_clock::now();
    cout << '.';
    duration<double> time_span_histo = time_span_solve_histo + duration_cast<duration<double>>(end - start);

    //////////////////////////////

    time_file << time_span_uds.count() << '\t'
              << time_span_wds.count() << '\t'
              << time_span_histo.count() << '\t'<< endl;

    ////////////////////////////// Estimate //////////////////////////////

    for (int param = 0; param < PARAM_UB; ++ param) {

      ////////////////////////////// UDS Estimate //////////////////////////////

      start = high_resolution_clock::now();
      RELATION_TYPE uds_sketch_filtered = filter_table(uds_sketch,param);
      unsigned int uds_sketch_distinct_count = get_distinct_count(uds_sketch_filtered);
      double uds_estimator = uds_sketch_distinct_count * pow(2, l_ds);
      uds_sketch_filtered.clear();
      end = high_resolution_clock::now();

      duration<double> time_span_uds = duration_cast<duration<double>>(end - start);

      ////////////////////////////// WDS Estimate //////////////////////////////

      start = high_resolution_clock::now();
      double wds_estimator = 0.0;
      for (auto it = wds_sketch.cbegin(); it != wds_sketch.cend(); ++it) {
        RELATION_TYPE Ra = it -> second;
        RELATION_TYPE Ra_filtered = filter_table(Ra,param);
        if (Ra_filtered.size() > 0) {
          wds_estimator += 1.0 / p_i[it -> first];
        }
        Ra_filtered.clear();
      }
      // wds_estimator = wds_estimator + (D - M) / 2.0;
      end = high_resolution_clock::now();

      duration<double> time_span_wds = duration_cast<duration<double>>(end - start);

      ////////////////////////////// HISTO Estimate //////////////////////////////

      start = high_resolution_clock::now();
      double histo_estimator = 0.0;
      for (auto it = histo_sketch.cbegin(); it != histo_sketch.cend(); ++it) {
        RELATION_TYPE Ra = it -> second;
        RELATION_TYPE Ra_filtered = filter_table(Ra,param);
        if (Ra_filtered.size() > 0) {
          histo_estimator += 1.0 / p_i_histo[it -> first];
        }
        Ra_filtered.clear();
      }
      histo_estimator = histo_estimator + (D - M_histo) / 2.0;
      end = high_resolution_clock::now();

      duration<double> time_span_histo = duration_cast<duration<double>>(end - start);


      ////////////////////////////// Exact //////////////////////////////

      start = high_resolution_clock::now();

      RELATION_TYPE filtered_data = filter_table(data,param);
      unsigned int filtered_size = filtered_data.size();
      unsigned int filtered_distinct_count = get_distinct_count(filtered_data);
      filtered_data.clear();
      end = high_resolution_clock::now();

      duration<double> time_span_exact = duration_cast<duration<double>>(end - start);

      time_file << time_span_uds.count() << '\t'
                << time_span_wds.count() << '\t'
                << time_span_histo.count() << '\t'
                << time_span_exact.count() << '\t' << endl;

      exact_file << filtered_distinct_count << '\t';
      ub_file << min(filtered_size, D) << '\t';
      uds_file << uds_estimator << '\t';
      wds_file << wds_estimator << '\t';
      histo_file << histo_estimator << '\t';
   }
    exact_file << endl;
    ub_file << endl;
    uds_file << endl;
    wds_file << endl;
    histo_file << endl;
    cout << endl;
    }
  exact_file.close();
  ub_file.close();
  uds_file.close();
  wds_file.close();
  histo_file.close();
  time_file.close();
  return 0;
}
