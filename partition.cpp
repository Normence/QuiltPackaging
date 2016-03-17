/*
   Copyright 2016 Peking University
   Author(s): Wenyu Luo <normence@gmail.com>
              Guojie Luo <gluo@pku.edu.cn>
*/

#include <iostream>
#include <fstream>
#include <string>
#include <boost/multi_array.hpp>
#include <boost/cstdlib.hpp>

using namespace std;
  
typedef boost::multi_array<double,5> DblArr5d;
typedef boost::multi_array<int,5> IntArr5d;
typedef boost::multi_array<int,4> IntArr4d;

/** (1 + Ac * d / alpha) ^ (-alpha) */
double area_yield(double area, double df_dnst, double alpha) {
  double yield = pow(1 + area*df_dnst/alpha, -alpha);
  return yield;
}

/** (1 - f_wire) ^ count */
double wire_yield(double f_wire, int count) {
  return pow(1 - f_wire, count);
}

/** find a recursive slicing partition with an optimal yield
  @opt    : the optimal yield
  @pos    : the cut position
  @w_size : #col of the partition grids
  @h_size : #row of the partition grids
  @act_w  : the width of the routing region
  @act_h  : the height of the routing region
 */
int Partition(DblArr5d& opt, IntArr5d& pos, IntArr4d& wire_count,
    int w_size, int h_size, int act_w, int act_h) {
  int L = 3; // maximum level of partition: an initial guess
  opt.resize(boost::extents[L+1][w_size][h_size][w_size+1][h_size+1]);
  pos.resize(boost::extents[L+1][w_size][h_size][w_size+1][h_size+1]);

  int x = 0, y = 0, w = 0, h = 0; // rectangle (x,y,w,h)
  int lx = 0;                     // lx: cut position at x=lx
  int ly = 0;                     // ly: cut position at y=ly
  int p = 0;                      // encoding of the cut position

  double unit;    // a unit length (in meters) of the routing data
  double df_dnst; // defect density
  double alpha;   // alpha (clustering parameter)
  double f_wire;  // failure rate of a cross-die wire

  cout << "Input the unit length (m): ";
  cin >> unit;
  cout << "Input the defect density: ";
  cin >> df_dnst;
  cout << "Input the clustering parameter: ";
  cin >> alpha;
  cout << "Input the wire-failure rate: ";
  cin >> f_wire;

  auto sqm = [unit,act_w,act_h,w_size,h_size](int grids) { // square meters
    double rsl_w = unit * act_w / double(w_size); // the width (m) of a partition grid
    double rsl_h = unit * act_h / double(h_size); // the height (m) of a partition grid
    return rsl_w * rsl_h * grids;
  };

  // traverse every possible rectangle
  for (h = 1; h <= h_size; ++h) {
    for (w = 1; w <= w_size; ++w) {
      for (y = 0; y <= h_size-h; ++y) {
        for (x = 0; x <= w_size-w; ++x) {
          double area = sqm(w*h);
          opt[0][x][y][w][h] = area_yield(area, df_dnst, alpha);
        }
      }
    }
  }

  int l = 1; // level of partition
  while (true) {
    while (true) {

      for (h = 1; h <= h_size; ++h) {
        for (w = 1; w <= w_size; ++w) {
          for (y = 0; y <= h_size-h; ++y) {
            for (x = 0; x <= w_size-w; ++x) {

              opt[l][x][y][w][h] = opt[l-1][x][y][w][h];
              pos[l][x][y][w][h] = -1; // no partition

              for (lx = x+1; lx <= x+w-1; ++lx) {
                p = lx - (x+1);
                double yield =
                  wire_yield(f_wire, wire_count[1][h][lx][y]) *
                  min(opt[l-1][x][y][lx-x][h], opt[l-1][lx][y][w-(lx-x)][h]);

                if (yield >= opt[l][x][y][w][h]) {
                  opt[l][x][y][w][h] = yield;
                  pos[l][x][y][w][h] = p;
                }
              }

              for (ly = y+1; ly <= y+h-1; ++ly) {
                p = (w-1) + ly - (y+1);  
                double yield =
                  wire_yield(f_wire, wire_count[0][w][x][ly]) *
                  min(opt[l-1][x][y][w][ly-y], opt[l-1][x][ly][w][h-(ly-y)]);

                if (yield >= opt[l][x][y][w][h]) {
                  opt[l][x][y][w][h] = yield;
                  pos[l][x][y][w][h] = p;
                }
              }
            }
          }
        }
      }

      cout << "level " << l << " :"
        << " yield = " << opt[l][0][0][w_size][h_size]
        << " pos = " << pos[l][0][0][w_size][h_size]
        << endl;

      if (l == L || pos[l][0][0][w_size][h_size] == -1)
        break;
      else
        l += 1;
    }
    if (pos[l][0][0][w_size][h_size] == -1)
      break;

    L += 3; // the next guess of the maximum level of partition
    opt.resize(boost::extents[L+1][w_size][h_size][w_size+1][h_size+1]);
    pos.resize(boost::extents[L+1][w_size][h_size][w_size+1][h_size+1]);
  }

  cout << "N x N = " << w_size << " x " << h_size << endl;
  cout << "defect density = " << df_dnst << "/m^2" << endl;
  cout << "clustering parameter = " << alpha << endl;
  cout << "wire failure rate = " << f_wire << endl;
  cout << "optimal yield = " << opt[l-1][0][0][w_size][h_size] << " at level " << l-1 << endl;

  return l-1;
}

/** read the graph file in the ISPD08 contest format
  @file_gr : the graph file
  @w       : the width of the routing region
  @h       : the height of the routing region
  @x       : the x coordinate of the lower-left corner
  @y       : the y coordinate of the lower-left corner
  */
bool ReadGraph(string file_gr, int &w, int &h, int &x, int &y) {
  ifstream fin(file_gr, ios::in);
  if ( ! fin.is_open()) return false;

  string no_use;
  fin >> no_use;

  int w1; // #col of the routing grids
  int h1; // #row of the routing grids
  fin >> w1 >> h1;

  fin.ignore(255, '\n');
  fin.ignore(255, '\n');
  fin.ignore(255, '\n');
  fin.ignore(255, '\n');
  fin.ignore(255, '\n');
  fin.ignore(255, '\n');

  fin >> x >> y;
  
  int w2; // the width of a routing grid
  int h2; // the height of a routing grid
  fin >> w2 >> h2;

  w = w1 * w2;
  h = h1 * h2;

  fin.close();

  return true;
}

/** is relatively prime */
bool is_rlt_prime(int a, int b) {
  int temp;
  if(a < b){
    temp = a;
    a = b;
    b = temp;
  }

  while((a %= b) != 0){
    temp = a;
    a = b;
    b = temp;
  }

  return (b == 1);
}

/** find a proper N' based on the given N */
int GetResolution(int w, int h, int N) {
  while (!is_rlt_prime(w, N) || !is_rlt_prime(h, N))
    N++;

  return N;
}

/** read the routing file in the ISPD08 contest format
  @file_rt    : the routing file
  @wire_count : the number of wires across edges in the partition grids
  @N          : #row & #col of the partition grids
  @act_w      : the width of the routing region
  @act_h      : the height of the routing region
  @ini_x      : the x coordinate of the lower-left corner
  @ini_y      : the y coordinate of the lower-left corner
  */
bool ReadRouting(string file_rt, IntArr4d &wire_count, int N,
                 int act_w, int act_h, int ini_x, int ini_y) {
  ifstream fin(file_rt, ios::in);
  if ( ! fin.is_open()) return false;

  double rsl_w = act_w / double(N); // the width of a partition grid
  double rsl_h = act_h / double(N); // the height of a partition grid

  int x1, y1, x2, y2;
  char ch;

  fin.get();

  while(!fin.eof()){
    fin.ignore(255, '\n');

    fin.get(ch);
    while (ch != '!') {
      fin >> x1;
      fin.get();
      fin >> y1;
      fin.ignore(5); // assume z in {0,1,...,9}
      fin >> x2;
      fin.get();
      fin >> y2;
      fin.ignore(255, '\n');

      // metric transformation
      //   before: (x,y) coordinates in the routing region
      //   after : col and row in the partition grids
      auto xform = [rsl_w,rsl_h,ini_x,ini_y](int& x, int& y) {
        x = int((x - ini_x) / rsl_w);
        y = int((y - ini_y) / rsl_h);
      };
      xform(x1, y1);
      xform(x2, y2);

      //if (x1 > x2 || y1 > y2) { swap(x1,x2); swap(y1,y2); }

      if (x1 < x2) {
        assert(y1 == y2);
        // horizontal wire, vertical edge
        for (int i=x1+1; i<=x2; i++)
          wire_count[1][1][i][y1]++;
      }
      else if (y1 < y2) {
        assert(x1 == x2);
        // vertical wire, horizontal edge
        for(int i=y1+1; i<=y2; i++)
          wire_count[0][1][x1][i]++;
      } else {
        assert(x1 == x2 && y1 == y2);
      }

      fin.get(ch);
    }

    fin.ignore(255, '\n');
    fin.get();
  }

  fin.close();

  for (int l=2; l<=N; l++) {
    // horizontal edges
    for (int x=0; x<=(N-l); x++)
      for (int y=1; y<=N-1; y++)
        wire_count[0][l][x][y] = wire_count[0][l-1][x][y] + wire_count[0][1][x+l-1][y];

    // vertical edges
    for(int y=0; y<=(N-l); y++)
      for(int x=1; x<=N-1; x++)
        wire_count[1][l][x][y] = wire_count[1][l-1][x][y] + wire_count[1][1][x][y+l-1];
  }

  return true;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: partition <graph> <routing>" << endl;
    return boost::exit_failure;
  }
  string file_gr = argv[1]; // the input (graph) file in the ISPD08 contest format
  string file_rt = argv[2]; // the output (routing) file in the ISPD08 contest format

  // parser begin //////////////////////////////////////////////////////////////
  int act_w; // the width of the routing region
  int act_h; // the height of the routing region
  int ini_x; // the x coordinate of the lower-left corner
  int ini_y; // the y coordinate of the lower-left corner
  cout << "Reading from the file '" << file_gr << "'..." << endl;
  if ( ! ReadGraph(file_gr, act_w, act_h, ini_x, ini_y))
    return boost::exit_failure;
  cout << "Succeed." << endl;

  int N; // #row & #col of the partition grids
  cout << "Input resolution N: ";
  cin >> N;
  N = GetResolution(act_w, act_h, N);
  cout << "The resolution N is adjusted to: " << N << endl;

  // the number of wires across an edge in the partition grids
  //   1st dim: 0 - horizontal edge; 1 - vertical edge
  //   2nd dim: the length of the edge in {1, 2, ..., N}
  //   3rd dim: the x coordinate of the left (bottom) end in {0, 1, ..., N}
  //   4th dim: the y coordinate of the left (bottom) end in {0, 1, ..., N}
  IntArr4d wire_count(boost::extents[2][N+1][N+1][N+1]);
  std::fill(wire_count.data(), wire_count.data() + wire_count.num_elements(), 0);
  
  cout << "Reading from the file '" << file_rt << "'..." << endl;
  // count the number of wires across any edge
  if ( ! ReadRouting(file_rt, wire_count, N, act_w, act_h, ini_x, ini_y))
    return -1;
  cout << "Succeed." << endl;
  // parser end ////////////////////////////////////////////////////////////////

  typedef boost::general_storage_order<5> storage;
  DblArr5d::size_type ordering_d5[] = {0, 4, 3, 2, 1};
  IntArr5d::size_type ordering_i5[] = {0, 4, 3, 2, 1};
  bool ascending[] = {true, true, true, true, true};

  // optimal yield (opt) and cut position (pos) of a given rectangle
  //   1st dim: number of levels for partition
  //   2nd dim: x coordinate of the lower-left corner in {0, 1, ..., N-1}
  //   3rd dim: y coordinate of the lower-left corner in {0, 1, ..., N-1}
  //   4th dim: the width of the rectangle in {1, 2, ..., N}
  //   5th dim: the height of the rectangle in {1, 2, ..., N}
  DblArr5d opt(
      boost::extents[2][N][N][N+1][N+1],
      storage(ordering_d5, ascending));
  IntArr5d pos(
      boost::extents[2][N][N][N+1][N+1],
      storage(ordering_i5, ascending));

  // the optimal level for partition
  int l = Partition(opt, pos, wire_count, N, N, act_w, act_h);
  
  return boost::exit_success;
}
