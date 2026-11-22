// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// items: 1-based item index per element (all sets concatenated)
// lens : length of each set
// n_items: number of distinct items
// method : 0 = overlap coefficient, 1 = Jaccard index
// [[Rcpp::export]]
NumericMatrix set_similarity_cpp(IntegerVector items, IntegerVector lens, int n_items, int method = 0, int nthreads = 1) {
  const int n = lens.size();
  std::vector<int> start(n + 1, 0);
  for (int i = 0; i < n; ++i) start[i + 1] = start[i] + lens[i];

  // 1) dedup within each set -> set -> items CSR
std::vector<int> mark(n_items, -1), sp(n + 1, 0), si;
si.reserve(items.size());
for (int i = 0; i < n; ++i) {
  for (int k = start[i]; k < start[i + 1]; ++k) {
    int it = items[k];
    if (it == NA_INTEGER) continue;
    --it;
    if (mark[it] != i) { mark[it] = i; si.push_back(it); }
  }
  sp[i + 1] = si.size();
}

// 2) item -> sets CSR (sets ascending)
std::vector<int> ip(n_items + 1, 0);
for (int it : si) ++ip[it + 1];
for (int t = 0; t < n_items; ++t) ip[t + 1] += ip[t];
std::vector<int> is(si.size()), pos(ip.begin(), ip.end() - 1);
for (int i = 0; i < n; ++i)
  for (int k = sp[i]; k < sp[i + 1]; ++k) is[pos[si[k]]++] = i;

NumericMatrix out(n, n);
double* o = out.begin();

// 3) for each set i: count shared items with sets j > i, then score
#pragma omp parallel num_threads(nthreads)
{
  std::vector<int> acc(n, 0), touched;
  touched.reserve(n);
  #pragma omp for schedule(dynamic, 16)
  for (int i = 0; i < n; ++i) {
    const int ni = sp[i + 1] - sp[i];
    if (ni == 0) continue;
    o[i + (size_t)i * n] = 1.0;
    for (int k = sp[i]; k < sp[i + 1]; ++k) {
      const int it = si[k];
      // sets are ascending: skip to those > i
      int a = ip[it], b = ip[it + 1];
      // binary search first index with is[] > i
      int lo = a, hi = b;
      while (lo < hi) { int m = (lo + hi) >> 1; if (is[m] <= i) lo = m + 1; else hi = m; }
      for (int q = lo; q < b; ++q) {
        const int j = is[q];
        if (acc[j]++ == 0) touched.push_back(j);
      }
    }
    for (int j : touched) {
      const int nj = sp[j + 1] - sp[j];
      const int c = acc[j];
      const double v = (method == 0)
      ? (double)c / (ni < nj ? ni : nj)   // |A∩B| / min(|A|,|B|)
      : (double)c / (ni + nj - c);        // |A∩B| / |A∪B|
        o[j + (size_t)i * n] = v;
        o[i + (size_t)j * n] = v;
        acc[j] = 0;
    }
    touched.clear();
  }
}
return out;
}
