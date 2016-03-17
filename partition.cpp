#include <iostream>
#include <fstream>
#include "boost/multi_array.hpp"
#include <string>

using namespace std;
	
typedef boost::multi_array<double, 5> fd_dbl;
typedef boost::multi_array<int, 5> fd_int;
typedef boost::multi_array<int, 4> fod_int;
typedef fd_dbl::index index_dbl; // not sure functions
typedef fd_int::index index_int;
typedef fod_int::index index_int;

extern double compute_base(int x, int y, int w, int h); // (1 + Ac * d / a) ^ (-a)
extern double compute_comp(); // (1 - f_{wire}) ^ N_{wire}



/**
 * N: resolution of partition
 */
int partition(fd_dbl &opt, fd_int &pos, int N, int x_size, int y_size, int w_size, int h_size) {

	int x = 0, y = 0, w = 0, h = 0; // rectangle (x,y,w,h)
	int l;							// iteration
	int lx = 0;                     // lx: cut position at x=lx
	int ly = 0;                     // ly: cut position at y=ly
	int p = 0;                      // encoding of the cut position
	double df_dnst, alpha, f_wire;	// defect density, alpha(clustering parameter), failure rate

	cout << "Input the defect density: ";
	cin >> df_dnst;
	cout << "Input alpha: ";
	cin >> alpha;
	cout << "Input the failure rate: ";
	cin >> f_wire;

	for (h = 1; h <= h_size; ++h) {  // traverse every possible rectangle
		for (w = 1; w <= w_size; ++w) {
			for (y = 0; y <= h_size-h; ++y) {
				for (x = 0; x <= w_size-w; ++x) {
					opt[x][y][w][h][0] = compute_base(x, y, w, h); // TO-DO
				}
			}
		}
	}

	int L = 3; // default initial value
	l = 1;

	while(1){
		for (; l <= L; ++l) {

			for (h = 1; h <= h_size; ++h) {
				for (w = 1; w <= w_size; ++w) {
					for (y = 0; y <= h_size-h; ++y) {
						for (x = 0; x <= w_size-w; ++x) {

							opt[x][y][w][h][l] = opt[x][y][w][h][l-1];
							pos[x][y][w][h][l] = -1; // no partition

							for (lx = x+1; lx <= x+w-1; ++lx) {
								p = lx - (x+1);
								double yield = compute_comp(); // TO-DO

								if (yield >= opt[x][y][w][h][l]) {
									opt[x][y][w][h][l] = yield;
									pos[x][y][w][h][l] = p;
								}
							}

							for (ly = y+1; ly <= y+h-1; ++ly) {
								p = (w-1) + ly - (y+1);  
								double yield = compute_comp(); // TO-DO

								if (yield >= opt[x][y][w][h][l]) {
									opt[x][y][w][h][l] = yield;
									pos[x][y][w][h][l] = p;
								}
							}
						}
					}
				}
			}

			if(pos[x][y][w][h][l] == -1)
				break;
		}

		L *= 2;
		opt.resize(extents[x_size][y_size][w_size+1][h_size+1][L]);
		pos.resize(extents[x_size][y_size][w_size+1][h_size+1][L]);
	}

	return l;

}

/** read the input file in the ISPD08 contest format
  @file_gr : the input file
  @w       : the width of the routing region
  @h       : the height of the routing region
  @x       : the x coordinate of the lower-left corner
  @y       : the y coordinate of the lower-left corner
  */
bool ReadGraph(string file_gr, int &w, int &h, int &x, int &y) {
	ifstream fin(file_gr, ios::in);
	if (!fin.is_open()) return false;

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

/** find a proper N' start from the given N */
int GetResolution(int w, int h, int N) {
	while (!is_rlt_prime(w, N) || !is_rlt_prime(h, N))
		N++;

	return N;
}

bool getData(int argc, char *argv, fod_int &com, int N, const double x_rsl, const double y_rsl, int ini_x, int ini_y)
{
	ifstream fin(argv[2], ios::in);
	if(!fin.is_open())
		return false;

	int x1, y1, x2, y2;
	char ch;

	fin.get();

	while(!fin.eof()){
		fin.ignore(255, '\n');

		fin.get(ch);
		while(ch != '!'){
			fin >> x1;
			fin.get();
			fin >> y1;
			fin.ignore(5);
			fin >> x2;
			fin.get();
			fin >> y2;
			fin.ignore(255, '\n');

			x1 -= ini_x;
			y1 -= ini_y;
			x2 -= ini_x;
			y2 -= ini_y;
			x1 = (int)(x1/x_rsl);  //metric transformation (from standard to resolution N)
			x2 = (int)(x2/x_rsl);
			y1 = (int)(y1/y_rsl);
			y2 = (int)(y2/y_rsl);

			if(x1 < x2){
				for(int i=x1+1; i<=x2; i++)
					com[1][1][i][y1]++;
			}
			else if(y1 < y2){
				for(int i=y1+1; i<=y2; i++)
					com[0][1][x1][i]++;
			}

			fin.get(ch);
		}

		fin.ignore(255, '\n');
		fin.get();
	}

	fin.close();

	for(int l=2; l<=N; l++) // horizontal
		for(int x=0; x<=(N-l); x++)
			for(int y=1; y<N; y++)
				com[0][l][x][y] = com[0][l-1][x][y] + com[0][1][x+l-1][y];

	for(int l=2; l<=N; l++) // vertical
		for(int y=0; y<=(N-l); y++)
			for(int x=1; x<N; x++)
				com[1][l][x][y] = com[1][l-1][x][y] + com[1][1][x][y+l-1];

	return true;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: partition <graph> <routing>" << endl;
    return EXIT_FAILURE;
  }
  string file_gr = argv[1]; // the input file in the ISPD08 contest format
  string file_rt = argv[2]; // the output file in the ISPD08 contest format

	// parser begin //////////////////////////////////////////////////////////////
	int act_w; // the width of the routing region
  int act_h; // the height of the routing region
  int ini_x; // the x coordinate of the lower-left corner
  int ini_y; // the y coordinate of the lower-left corner
	cout << "Reading from the file '" << file_gr << "'..." << endl;
	if ( ! ReadGraph(file_gr, act_w, act_h, ini_x, ini_y))
		return EXIT_FAILURE;
	cout << "Succeed." << endl;

	int N; // #row & #col of the partition grids
	cout << "Input resolution N: ";
  cin >> N;
  cout << "The resolution N is adjusted from " << N
      << " to " << (N = GetResolution(act_w, act_h, N))
      << endl;

	const double x_rsl = act_w / N; // the width of a partition grid
  const double y_rsl = act_h / N; // the height of a partition grid
	int x_size = N;
  int y_size = N;
  int w_size = N;
  int h_size = N;

  // the number of wires across an edge in the partition grids
  //   1st dim: 0 - horizontal edge; 1 - vertical edge
  //   2nd dim: the length of the edge in {1, 2, ..., N-1}
  //   3rd dim: the x coordinate of the left (bottom) end
  //   4th dim: the y coordinate of the left (bottom) end
	fod_int com(boost::extents[2][N][w_size+1][h_size+1]);
  std::fill(com.data(), com.data() + com.num_elements(), 0);
	
	cout << "Reading from the file '" << file_rt << "'..." << endl;
  // count the number of wires across any edge
	if ( ! getData(argc, argv, com, N, x_rsl, y_rsl, ini_x, ini_y))
		return -1;
	cout << "Succeed." << endl;
	// parser end ////////////////////////////////////////////////////////////////


	fd_dbl opt(boost::extents[x_size][y_size][w_size+1][h_size+1][L]); // optimal yield
	fd_int pos(boost::extents[x_size][y_size][w_size+1][h_size+1][L]); // cut position

	int l;  // best recursive level
	l = partition(opt, pos, N, x_size, y_size, w_size, h_size);
	
	return EXIT_SUCCESS;
}
