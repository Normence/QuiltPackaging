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
extern double compute_comp(); // (1 - f_{fire}) ^ N_{wire}


/**
 * N: resolution of partition
 * L: maximum recursive levels
 */
int partition(fd_dbl &opt, fd_int &pos,int N) {

	int x = 0, y = 0, w = 0, h = 0; // rectangle (x,y,w,h)
	int l;  // iteration
	int lx = 0;                     // lx: cut position at x=lx
	int ly = 0;                     // ly: cut position at y=ly
	int p = 0;                      // encoding of the cut position


	for (h = 1; h <= h_size; ++h) {  // traverse every possible rectangle
		for (w = 1; w <= w_size; ++w) {
			for (y = 0; y <= h_size-h; ++y) {
				for (x = 0; x <= w_size-w; ++x) {
					opt[x][y][w][h][0] = compute_base(x, y, w, h); // TO-DO
				}
			}
		}
	}

	for (l = 1; l < L; ++l) {

		for (h = 1; h <= h_size; ++h) {
			for (w = 1; w <= w_size; ++w) {
				for (y = 0; y <= h_size-h; ++y) {
					for (x = 0; x <= w_size-w; ++x) {

						opt[x][y][w][h][l] = opt[x][y][w][h][l-1];
						pos[x][y][w][h][l] = -1; // no partition
					
						for (lx = x+1; lx <= x+w-1; ++lx) {
							p   = lx - (x+1);
							double yield = compute_comp(); // TO-DO

							if (yield >= opt[x][y][w][h][l]) {
								opt[x][y][w][h][l] = yield;
								pos[x][y][w][h][l] = p;
							}
						}
						
						for (ly = y+1; ly <= y+h-1; ++ly) {
							p   = (w-1) + ly - (y+1);  
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

	return l;

}

bool getData(int argc, char *argv[], &w, &h, &x, &y) // get left-bottom (x, y), width and height
{
	if(argc != 3)
		return false;

	ifstream fin(argv[1], ios::in);
	if(!fin.is_open())
		return false;

	string no_use;
	fin >> no_use;

	int w1, h1;
	fin >> w1 >> h1;

	fin.ignore(255, '\n');
	fin.ignore(255, '\n');
	fin.ignore(255, '\n');
	fin.ignore(255, '\n');
	fin.ignore(255, '\n');
	fin.ignore(255, '\n');

	fin >> x >> y;
	
	int w2, h2;
	fin >> w2 >> h2;

	w = w1 * w2;
	h = h1 * h2;

	fin.close();

	return true;
}

bool is_rlt_prime(int a, int b)  // is relatively prime
{
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

	if(b == 1)
		return true;
	else
		return false;
}

int getResolution(int w, int h, int N)
{
	while(!is_rlt_prime(w, N) || !is_rlt_prime(h, N))
		N++;

	return N;
}

bool getData(int argc, char *argv, fod_int com, int N, double x_rsl, double y_rsl, int ini_x, int ini_y)
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

	return true;
}

int main(int argc, char *argv[])
{
	int act_w, act_h, ini_x, ini_y;
	if(!getData(argc, argv, act_w, act_h, ini_x, ini_y))
		return -1;

	int N;
	cout << "Input resolution N: ";
	cin >> N;
	N = getResolution(act_w, act_h, N);

	const double x_rsl = act_w / N, y_rsl = act_h / N;

	int x_size = N, y_size = N, w_size = N, h_size = N;

	fod_int com(boost::extents[2][N][w_size+1][h_size]+1); // [vertical/horizontal][length][x][y]
	for(int i=0; i<2; i++)
		for(int j=0; j<N; j++)
			for(int q=0; q<w_size+1; q++)
				for(int p=0; p<h_size+1; p++)
					com[i][j][q][p] = 0;

	if(!getData(argc, argv, com, N, x_rsl, y_rsl, ini_x, ini_y))
		return -1;

	// L...to-do
	fd_dbl opt(boost::extents[x_size][y_size][w_size+1][h_size+1][L]); // optimal yield
	fd_int pos(boost::extents[x_size][y_size][w_size+1][h_size+1][L]);  // cut position

	int l;  // best recursive level
	l = partition(opt, pos, N);
	
	return 0;
}
