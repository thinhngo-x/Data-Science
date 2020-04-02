#include <iostream>  // for cout
#include <fstream>  // for ifstream
#include <cfloat>  // for DBL_MAX
#include <cmath>  // for sqrt
#include <cstdlib>  // for rand, srand
#include <ctime>  // for rand seed
#include <cassert>  // for assertions
#include <algorithm>  // for sort
#include <ctime>  // for clock

const bool debug = false;  // debug flag, r
// remember to turn it back to false before submitting your code for automatic validation.

typedef double* point;  // point = array of coordinates (doubles)

typedef struct com{
  int c;
  com(int c){this->c = c;};
  bool operator ()(point a, point b){
    return a[c] < b[c];
  }
}com;

void print(point p, int dim) {
  std::cout << p[0];
  for ( int j = 1; j < dim; j++)
    std::cout << " " << p[j];
  std::cout << "\n";
}

void swap (point* P, int i, int j) {  // swap 2 points in data set
  point temp = P[i];
  P[i] = P[j];
  P[j] = temp;
}

double dist (point p, point q, int dim) {
  double res = 0;
  for(int i = 0; i < dim; i++)
    res += (p[i] - q[i]) * (p[i] - q[i]);
  res = std::sqrt(res);
  return res;
}

void test_dist(){
  int dim = 4;
  point P[2];
  P[0] = new double[dim];
  P[1] = new double[dim];
  for(int i = 0; i < dim; i++){
    P[0][i] = 0.0d;
    P[1][i] = 1.0d;
  }
  std::cout<<dist(P[0], P[1], dim)<<std::endl;
}

int linear_scan (point q, int dim, point* P, int n) {
  double mindist = DBL_MAX;
  int res = -1;
  for(int i = 0; i < n; i++){
    double tmp = dist(q, P[i], dim);
    if(tmp < mindist){
      res = i;
      mindist = tmp;
    }
  }
  return res;
}

void test_linear_scan(){
  int dim = 10;
  int n = 4;
  point P[4];
  for(int i = 0; i < n; i++)
    P[i] = new double[dim];
  for(int j = 0; j < dim; j++)
    for(int i = 0; i < n; i++){
      P[i][j] = 0 + i;
    }
  std::cout<<linear_scan(P[0], dim, P, n)<<std::endl;
}

double computeMedian(point* P, int start, int end, int c) {
  double cth_value[end-start];
  for(int i = 0; i < end-start; i++){
    cth_value[i] = P[i+start][c];
  }
  std::sort(cth_value, cth_value + end-start);
  // std::cout<<(int) ((end-start)/2)<<std::endl;
  return cth_value[(end-start)/2];
}

void test_computeMedian(){
  int dim = 10;
  int n = 4;
  point P[4];
  for(int i = 0; i < n; i++)
    P[i] = new double[dim];
  for(int j = 0; j < dim; j++)
    for(int i = 0; i < n; i++){
      P[i][j] = i+j;
    }
  std::cout<<computeMedian(P, 0, n, 0)<<std::endl;
}

int partition(point* P, int start, int end, int c, int dim) {
  // int i = start;
  // int j = end-1;
  // int p = (start+end)/2;
  // double med = computeMedian(P, start, end, c);
  // while(i <= j){
  //   if(P[i][c] == med && P[j][c] == med)
  //     break;
  //   if(P[i][c] >= med && P[j][c] <= med)
  //     swap(P, i, j);
  //   else if(P[i][c] < med)
  //     i++;
  //   else if(P[j][c] > med)
  //     j--;
  // }

  std::sort(P + start, P + end, com(c));
  int p = (start+end)/2;

  // std::sort(P+start, P+end, [](point* a, point* b){return a[c] < b[c]});
  // int p = (start+end)/2;
  if(debug){
    std::cout<<"After partition accord c="<<c<<" from start="<<start<<" to end="<<end<<std::endl;
    for(int i = start; i < end; i++)
      print(P[i], dim);
  }
  return p;
}

void test_partition(){
  int dim = 10;
  int n = 4;
  point P[4];
  for(int i = 0; i < n; i++)
    P[i] = new double[dim];
  for(int j = 0; j < dim; j++)
    for(int i = 0; i < n; i++){
      P[i][j] = n-i;
    }
  std::cout<<partition(P, 0, n, 0, dim)<<std::endl;
  for(int i = 0; i < n; i++){
    print(P[i], dim);
  }
}

typedef struct node {  // node for the kd-tree
  int c;
  double m;
  int p;
  node* left, *right;
} node;

node* create_node (int _p) {  // creates a leaf node
  node* init = new node;
  init->p = _p;
  init->left = NULL;
  init->right = NULL;
  if(debug){
    std::cout<<"Leaf Node: "<<init->p<<std::endl;
  }
  return init;
}

node* create_node (int _c, double _m, int _p, node* _left, node* _right) {  // creates an internal node
  node* init = new node;
  init->p = _p;
  init->c = _c;
  init->m = _m;
  init->left = _left;
  init->right = _right;
  if(debug){
    std::cout<<"Internal Node: "<<init->p<<std::endl;
  }
  return init;
}

node* build (point* P, int start, int end, int c, int dim) {
  // builds tree for sub-cloud P[start -> end-1]
  assert (end-start >= 0);
  if (debug)
    std::cout << "start=" << start << ", end=" << end << ", c=" << c << std::endl;
  if (end-start == 0)  // no data point left to process
    return NULL;
  else if (end-start == 1)  // leaf node
    return create_node (start);
  else {  // internal node
    if (debug) {
      std::cout << "array:\n";
      for (int i=start; i<end; i++)
        print(P[i],dim);
	// std::cout << P[i] << ((i==end-1)?"\n":" ");
    }
    // compute partition
    // rearrange subarray (less-than-median first, more-than-median last)
    int p = partition (P, start, end, c, dim);
    double m = P[p][c];
    // prepare for recursive calls
    int cc = (c+1)%dim;  // next coordinate
    return create_node (c, m, p,
		     build (P, start, p, cc, dim),
		     build (P, p+1, end, cc, dim));
  }
}

void test_build(){
  int n = 5;  // n points in R^{dim}
  int dim = 3;
  // // random input data points (uniformly sampled in unit cube)
  // srand (time(NULL));
  // point P[n];
  // for (int i=0; i<n; i++) {
  //   P[i] = new double [dim];
  //   for (int j=0; j<dim; j++)
  //     P[i][j] = (double)rand() / RAND_MAX;
  // }
  point P[n];
  for (int i=0; i<n; i++)
    P[i] = new double [dim];
  P[0][0] = 1.0d; P[0][1] = 7.0d; P[0][2] = 1.0d;
  P[1][0] = 2.0d; P[1][1] = 2.0d; P[1][2] = 2.0d;
  P[2][0] = 2.0d; P[2][1] = 1.0d; P[2][2] = 3.0d;
  P[3][0] = 2.0d; P[3][1] = 3.0d; P[3][2] = 1.0d;
  P[4][0] = 1.0d; P[4][1] = 2.0d; P[4][2] = 3.0d;
  for (int i=0; i<n; i++)
    print(P[i], dim);
  node* root = build(P, 0, n, 0, dim);
  for (int i=0; i<n; i++)
    print(P[i], dim);
}
  
void dsearch (node* n, point q, int dim, point* P, double& res, int& nnp) {
  if (n != NULL) {
    double tmp = dist(q, P[n->p], dim);
    if (res > tmp){
      nnp = n->p;  
      res = tmp;
    }
    if (n->left != NULL || n->right != NULL)  // internal node
      dsearch ( (q[n->c] <= n->m)?n->left:n->right, q, dim, P, res, nnp);
  }
}

void bsearch (node* n, point q, int dim, point* P, double& res, int& nnp) {
  if(n != NULL){
    double tmp = dist(q, P[n->p], dim);
    if(res > tmp){
      nnp = n->p;
      res = tmp;
    }
    if(q[n->c] + res > n->m && n->right != NULL)
      bsearch(n->right, q, dim, P, res, nnp);
    if(q[n->c] - res <= n->m && n->left != NULL)
      bsearch(n->left, q, dim, P, res, nnp);
  }
}



int main () {
  // test_dist();
  // test_linear_scan();
  // test_computeMedian();
  // test_partition();
  // test_build();

  const int dim = 4;  // dimension (hard-coded)
  int n = 20000;  // upper bound on number of data points in R^{dim}
  int nt = 1000;  // nt query points

  point P[n];
  char names[n][255];
  
  srand (time(NULL));
  
  for (int i=0; i<n; i++)
    P[i] = new double[dim];
  std::ifstream is("iris2.data");
  assert(is.is_open());
  for (n=0; is >> P[n][0]; n++) {
    for (int i=1; i<dim; i++)
      is >> P[n][i];
    is >> names[n];
  }
  std::cout << n << " observations" << std::endl;


/********* section 2 ***************/
  
  // //Some inputs / outputs
  // while (true){
  //   std::cout << "\nEnter your own measurements (4 space-separated real numbers in [0,10]):\n(type CTRL-D to exit): " << std::endl;
  //   if (std::cin.peek() == EOF) break;
  //   point query = new double [dim];
  //   std::cin >> query[0] >> query[1] >> query[2] >> query[3];
  //   double distb = DBL_MAX;
  //   int nn;
  //   nn = linear_scan (query, dim, P, n);
  //   // bsearch(tree, query, dim, P, distb, nn);
  //   std::cout << "\n  -> Your iris is of type " << names[nn] << std::endl;
  // }

  
/******* section 3 ********************************************/    

  // build kd-tree (warning: rearranges points in data set)
  std::cout << "building kd-tree..." << std::flush;
  node* tree = build (P, 0, n, 0, dim);
  std::cout << " done" << std::endl;
  

  double* d = new double[nt];  // recorded distances
  point* q = new point[nt];  // query points
  int* nnp = new int[nt];  // nearest neighbors

  // random query points
  for (int i=0; i<nt; i++) {
    q[i] = new double [dim];
    for (int j=0; j<dim; j++)
      q[i][j] = 10*(double)rand() / RAND_MAX;
  }


  // benchmarking, start and end time for chronos
  std::clock_t start, end;


  // linear scan
  std::cout << "Benchmarking linear scan..." << std::flush;
  start = std::clock();  
  for (int i=0; i<nt; i++) {
    int nn = linear_scan (q[i], dim, P, n);
    d[i] = dist(q[i], P[nn], dim);
    nnp[i] = nn;
    // std::cout << "nearest neighbor (linear scan) = P[" << nn
    // 	      << "] / distance = " << d[i] << std::endl; 
  }
  end = std::clock();  
  std::cout << " done (avg time = "
	    << (end - start)/nt
	    << " us)" << std::endl;


  // defeatist search
  std::cout << "Benchmarking defeatist search..." << std::flush;
  int nf = 0;
  start = std::clock();  
  for (int i=0; i<nt; i++) {
    double distd = DBL_MAX;
    int nn;
    dsearch(tree, q[i], dim, P, distd, nn);
    assert (distd >= d[i]);
    if (nn != nnp[i]) nf++;
    // std::cout << "nearest neighbor (defeatist search): " << nnp
    // 	      << " | distance = " << distd << std::endl; 
  }
  end = std::clock();  
  std::cout << " done (avg time = "
	    << (end - start)/nt
	    << " us, failure rate = " << (nf*100.0/nt)
	    << "\%)" << std::endl;
  

  // backtracking search
  std::cout << "Benchmarking backtracking search..." << std::flush;
  nf = 0;
  start = std::clock();  
  for (int i=0; i<nt; i++) {
    // std::cerr << i << " ";
    double distb = DBL_MAX;
    int nn;
    bsearch(tree, q[i], dim, P, distb, nn);
    assert (nn == nnp[i]); 
    if (distb != d[i]) nf++;
  }
  end = std::clock();  
  std::cout << " done (avg time = "
	    << (end - start)/nt
	    << " us, failure rate = " << (nf*100.0/nt)
	    << "\%)" << std::endl;
  
  return 0;
}
