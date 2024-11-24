#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <algorithm>

#include <sstream>

// #include <integration.h>
// #include <alglibinternal.h>
// #include <ap.h>

// using namespace alglib_impl;
using namespace std;

int main(int argc, char const *argv[])
{
	ifstream infile;

	int n_regions = 0;
	double tol = 0, alpha, max_iters;
	// alglib::ae_int_t n_mu;
	int n_mu; // temporary
	// define vectors for dynamic memory allocation
	std::vector<double> Lengths;
	std::vector<int> n_cells;
	std::vector<double> Sig_ts;
	std::vector<double> Sig_ss;
	std::vector<double> sources; // store sources
	int div_val;

	infile.open("input_file.txt");

	// attempt to open input file
	if (infile.is_open())
	{
		// assign variables here

		string line; // holds the current line in the input file

		while (getline(infile, line))
		{
			if (!(line[0] == '/' && line[1] == '/'))
			{
				if (tol == 0)
				{

					std::istringstream iss(line);			  // Create a string stream from the line
					iss >> tol >> max_iters >> n_mu >> alpha; // Extract the integer value
					continue;
				}
				else
				{
					std::istringstream iss(line); // Create a string stream from the line
					iss >> n_regions;			  // Extract the integer value
					break;
				}
			}
		}

		cout << "Computer settings: " << endl
			 << "	tol = " << tol << endl
			 << "	max_iters = " << max_iters << endl
			 << "	quadrature nodes = " << n_mu << endl
			 << "	alpha = " << alpha << endl
			 << endl;

		int i = 0; // initialize a counter

		while (getline(infile, line)) // test for a comment, denoted with '//'
		{

			// Skip comment lines
			if (line.empty() || (line[0] == '/' && line[1] == '/'))
			{
				continue;
			}

			istringstream iss(line); // convert line to string stream
			string token;			 // extract values separated by spaces

			while (iss >> token)
			{
				div_val = i/n_regions;
				cout << "div_val = " << div_val <<endl;
				if (div_val == 0)
				{
					Lengths.push_back(stod(token)); // 'stod' converts string to double
					cout << "DEBUG: Length = " << stod(token) << endl << endl;
				}
				else if (div_val == 1)
				{
					n_cells.push_back(stoi(token)); // 'stoi' converts string to int
				}
				else if (div_val == 2)
				{
					Sig_ts.push_back(stod(token));
				}
				else if (div_val == 3)
				{
					Sig_ss.push_back(stod(token));
				}
				else if (div_val == 4)
				{
					sources.push_back(stod(token));
				}

				i += 1;
				cout << "i = " << i << endl;
			}
		}

		infile.close(); // close the file to save resources
	}
	else
	{
		cout << "INPUT ERROR! FILE NOT FOUND!!!";
		return 1;
	}

	// Check for successful parsing
	cout << "Run settings: " << endl;
	cout << "Number of regions = " << n_regions << endl;
	for (int j = 0; j < n_regions; j++)
	{
		cout << "	Region " << j + 1 << ": " << endl
			 << "		Length = " << Lengths[j] << endl
			 << "		Cells  = " << n_cells[j] << endl
			 << "		Sig_t  = " << Sig_ts[j] << endl
			 << "		Sig_s  = " << Sig_ss[j] << endl
			 << "		Source = " << sources[j] << endl;
	}

	return 0;
}