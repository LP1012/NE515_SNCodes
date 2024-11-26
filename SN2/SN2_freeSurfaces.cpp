// To run code:
// g++ -I./src SN2_freeSurfaces.cpp ./src/*.cpp -o SN2_freeSurfaces && ./SN2_freeSurfaces

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <cmath>

#include <sstream>

#include <integration.h>
#include <alglibinternal.h>
#include <ap.h>

using namespace alglib_impl;
using namespace std;

double phi_forward(double Sigt, double dz, double mu_n, double lh_flux, double q)
{
    double den = 1 + Sigt * dz / (2 * mu_n);
    double num = lh_flux + dz * q / (2 * mu_n);
    return num / den;
}

double forward_flux(double alpha, double phi_n, double lh_flux)
{
    return 1 / (1 - alpha) * (phi_n - alpha * lh_flux);
}

double phi_backward(double Sigt, double dz, double mu_n, double rh_flux, double q)
{
    double den = 1 + Sigt * dz / (2 * abs(mu_n));
    double num = rh_flux + dz * q / (2 * abs(mu_n));
    return num / den;
}

double backward_flux(double alpha, double phi_n, double rh_flux)
{
    return 1 / (1 - alpha) * (phi_n - alpha * rh_flux);
}

double GQ_integrate(const std::vector<double> &f, const std::vector<double> &weights)
// Function evaluates the Gaussian quadrature using the function values and weights
{
    double sum = 0;
    for (int k = 0; k < f.size(); k++)
    {
        sum += f[k] * weights[k];
    }
    return sum;
}

// ----------------------------------------------------------------------
// Boundary Condition Functions -- ANALYTICAL
// ** These can be changed to whatever function is needed
// ----------------------------------------------------------------------
double Gamma_L(double mu, double mu_star)

{
    if (mu == mu_star)
    {
        return 1.0 / (2 * M_PI * mu_star);
    }
    else
    {
        return 0;
    }
}

double Gamma_R()
{
    return 0;
}

double BC_scaling(const std::vector<double> &mus, const std::vector<double> &weights, const std::vector<double> &gammas, string flag)
// Function calculates the scaling factor for the boundary fluxes
// Variable "flag" specifies whether the flux is on the right or left boundary; takes values of "R" and "L", respectively
{

    if (flag == "L")
    {

        int N = mus.size();
        int n = mus.size() / 2;
        double sum = 0;
        std::vector<double> fs;
        std::vector<double> new_weights;

        for (int k = n; k < N; k++)
        {
            fs.push_back(mus[k] * gammas[k]);
            new_weights.push_back(weights[k]);
        }

        double term = GQ_integrate(fs, new_weights);
        return 2 * M_PI * term;
    }
    else if (flag == "R")
    {
        int N = mus.size();
        int n = mus.size() / 2;
        double sum = 0;
        std::vector<double> fs;
        std::vector<double> new_weights;

        for (int k = n; k < N; k++)
        {
            fs.push_back(mus[k] * gammas[k]);
            new_weights.push_back(weights[k]);
        }

        double term = GQ_integrate(fs, new_weights);
        return 2 * M_PI * term;
    }
    else
    {
        printf("\nINVALID FLAG FOR SCALING FACTOR!!!\n");
        return -1;
    }
}

int find_region_index(int n, std::vector<int> num_cells)
// Function recursively calls itself to find the current region the index, n, is in
{
    // first compute the sum
    int sum = 0;
    for (int i = 0; i < num_cells.size(); i++)
    {
        sum += num_cells[i];
    }

    if (n < sum && num_cells.size() > 1)
    {
        num_cells.pop_back();
        return find_region_index(n, num_cells);
    }
    else if (n >= sum)
    {
        return num_cells.size();
    }
    else if (num_cells.size() == 1)
    {
        return 0;
    }

    cout << "ERROR in find_region_index --> No value returned" << endl;
    cout << "n = " << n << endl;
    return -1;
}

// ----------------------------------------------------------------------

int main(int argc, char const *argv[])
{
    printf("\n****************************************************************\n");
    printf("                 SN2 Free Surface CODE EXECUTION\n");
    printf("****************************************************************\n\n");
    std::vector<double> fmflux;      // allocate memory for storing cell fluxes -- "forward marching flux"
    std::vector<double> bmflux;      // allocate memory for storing cell fluxes -- "backward marching flux"
    int region_index;                // hold index number of current iteration -- to be used later
    double dz, Sig_t, Sig_s, source; // allocate memory for holding a temporary variable for each regions -- to be used in looops
    double running_error;
    std::vector<double> error_vec;
    string output_filename;

    bool zero_dirichlet_bool = true; // as the name of the program suggests, this flag specifies zero boundary conditions on slab ends

    // BEGIN INPUT TEXT FILE ---------------------------------------------------------------------------------------

    printf("Input data filename (no extension): ");
    cout << endl;

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
    int div_val;

    string infile_name;
    cin >> infile_name;

    string infile_with_extension;
    infile_with_extension = infile_name + ".txt";

    infile.open(infile_with_extension);

    printf("Importing data from text file: ");

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

                    std::istringstream iss(line);             // Create a string stream from the line
                    iss >> tol >> max_iters >> n_mu >> alpha; // Extract the integer value
                    continue;
                }
                else
                {
                    std::istringstream iss(line); // Create a string stream from the line
                    iss >> n_regions;             // Extract the integer value
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
            string token;            // extract values separated by spaces

            while (iss >> token)
            {
                div_val = i / n_regions;
                if (div_val == 0)
                {
                    Lengths.push_back(stod(token)); // 'stod' converts string to double
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

                i += 1;
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
             << "		Sig_s  = " << Sig_ss[j] << endl;
    }
    printf("\nData import complete.\n");

    // END INPUT TEXT FILE --------------------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // Calculate total number of cells and step sizes
    // ----------------------------------------------------------------------
    printf("Calculating total number of cells... ");
    int tot_n_cells = 0;     // store total number of cells
    std::vector<double> dzs; // store step size
    for (int i = 0; i < n_regions; i++)
    {
        tot_n_cells += n_cells[i];
        dzs.push_back((double)Lengths[i] / (double)n_cells[i]);
    }
    printf("Done.\n"
           "   Total number of cells = %d\n\n",
           tot_n_cells);

    // ----------------------------------------------------------------------
    // Gauss Quadrature data
    // ----------------------------------------------------------------------

    // Variables to store the output
    alglib::ae_int_t info;          // To hold the status code
    alglib::real_1d_array mus_init; // To store nodes
    alglib::real_1d_array w;        // To store weights
    std::vector<double> w_vec;      // Store weights as vector
    double mus[n_mu] = {};          // Store discrete ordinates in array
    std::vector<double> mu_vec;     // Store DO's in vector (for integration)

    // Call the Gauss-Legendre quadrature generator
    printf("Calculating GQ weights and discrete ordinates... ");
    alglib::gqgenerategausslegendre(n_mu, info, mus_init, w);
    printf("Done.\n"
           "Casting GQ data to a vector (weights) and array/vector (scattering angle)...");

    // Cast weights into a vector and mu's into an array and a vector
    for (int i = 0; i < n_mu; i++)
    {
        w_vec.push_back(w[i]);
        mus[i] = mus_init[i];
        mu_vec.push_back(mus_init[i]);
    }
    printf("Done.\n"
           "Initializing variables and iterating over cell sizes...\n");
    std::cout << "\nGQ Weights:  w = " << w.tostring(4);
    std::cout << "\nGQ Nodes:   mu = " << mus_init.tostring(4) << "\n";

    double mu_n, z_i, front_coeff;        // Allocate memory for variables in loops
    double fFlux, bFlux;                  // Corresponding to "front flux" and "back flux", respectively
    std::vector<double> angular_flux_col; // Vector to store column of angular flux values (corresponding to the same point in space)

    // ----------------------------------------------------------------------
    // Define boundary conditions
    // ----------------------------------------------------------------------

    printf("\nAssigning boundary conditions...\n");
    // Define analytical fluxes
    std::vector<double> gammas_L, gammas_R;

    double mu_star = mus[n_mu - 1]; // Choose the largest mu value just for this application
    for (int k = 0; k < n_mu; k++)
    {
        gammas_L.push_back(Gamma_L(mus[k], mu_star));
        gammas_R.push_back(Gamma_R());
    }

    printf("    Analytical fluxes set.\n");

    if (zero_dirichlet_bool)
    {
        printf("    Vacuum boundary conditions specified. Scaling factor calculation ignored.\n"
               "    Assigning zero fluxes on boundary... ");

        for (int k = 0; k < n_mu; k++)
        {
            gammas_L[k] = 0;
            gammas_R[k] = 0;
        }
        printf("Done.\n");
    }
    else
    {
        printf("    Calculating scaling factor... ");
        // Demand current to be 1
        double BC_scaling_factor_left = BC_scaling(mu_vec, w_vec, gammas_L, "L"); // calculated LH scaling factor
        printf("Done.\n"
               "    LH Scaling factor = %.4f\n",
               BC_scaling_factor_left);

        // *** UNCOMMENT THIS FOR NONZERO RH BOUNDARY CONDITION ***
        // double BC_scaling_factor_right = BC_scaling(mu_vec,w_vec,gammas_R,"R"); // calculate RH scaling factor

        printf("    Scaling numerical flux values... ");
        for (int k = 0; k < n_mu; k++)
        {
            gammas_L[k] *= (1 / BC_scaling_factor_left);

            // *** UNCOMMENT THIS FOR NONZERO RH BOUNDARY CONDITION ***
            // gammas_R[k] *= (1/BC_scaling_factor_right);
        }
        printf("Done.\n");
    }

    // ----------------------------------------------------------------------
    // Begin looping over number of cells in spatial coordinates
    // ----------------------------------------------------------------------

    printf("\n---------------------------------------------\n"
           "Begin looping over cells...\n"
           "---------------------------------------------\n\n");

    // Create 2D arrays to hold the scalar flux at different mu and z values
    // mu is varied in rows, while z is varied by column
    printf("    Initializing angular flux arrays... ");
    double psi_current[n_mu][tot_n_cells] = {};
    double psi_old[n_mu][tot_n_cells] = {};
    printf("Done.\n");

    // Must initialize q because it is used in the first iteration.
    // When compared to the notes, this will be q^(n-1), and code will simply make a direct update
    // instead of creating a new vector.
    printf("    Intializing scattering source terms... ");
    // Because of angular dependence, both must be 2D arrays.
    // Same structure as the angular flux (psi) values
    double q[n_mu] = {};
    double S[n_mu] = {};

    for (int i = 0; i < n_mu; i++)
    // loop over angles
    {
        mu_n = mus[i];
        S[i] = 1 / (2.0 * M_PI * w_vec[i] * (double)n_mu); // pull the source value based on the region (assumes constant source in a region)
        q[i] = 0;                                          // for the first iteration, we assume a 0 scalar flux everywhere
    }
    printf("Done.\n");

    // check source
    double s_sum = 0;
    for (int i = 0; i < n_mu; i++)
    {
        s_sum += 2 * M_PI * w_vec[i] * S[i];
    }
    cout << "    Source normalization check: Source integral = " << s_sum <<endl;
    

    printf("    Initialize computational error... ");
    double error = 1; // Dummy value to start with
    printf("Done.\n");
    int iteration_num = 0; // Keep track of the number of iterations passed

    // Cast scalar fluxes as vectors to keep with the data type in the functions written
    std::vector<double> scalar_flux_current;
    std::vector<double> scalar_flux_old;

    printf("    Beginning source iteration...\n");

    while (error > tol)
    {
        // ----------------------------------------------------------------------
        // Iterate over scattering angles (discrete ordinates)
        // ----------------------------------------------------------------------
        for (int i = n_mu - 1; i > -1; i--) // Step from max(mu_n) to min(mu_n)
        {
            // Update mu_n value
            mu_n = mus[i];

            if (mu_n > 0)
            {
                bmflux.clear();
                bmflux.push_back(gammas_L[i]);
                for (int j = 0; j < tot_n_cells; j++)
                {
                    region_index = find_region_index(j, n_cells); // determines the spatial region we are in
                    dz = dzs[region_index];
                    Sig_t = Sig_ts[region_index];
                    // Forward sweep in z
                    bFlux = bmflux.back();
                    front_coeff = 1 / (1 + Sig_t * dz / (2 * mu_n));
                    psi_current[i][j] = front_coeff * (bFlux + dz * q[i] / (2 * mu_n));

                    // This handles when the flux goes negative.
                    // Take the incoming flux to be equal to the previous cell-center flux (alpha=0)
                    if (psi_current[i][j] < 0)
                    {
                        bmflux.erase(std::find(bmflux.begin(), bmflux.end(), bFlux));          // Delete the final flux in bmflux, the value that made us go negative
                        bFlux = psi_current[i][j - 1];                                         // reassign flux for alpha=0 only in this cell
                        bmflux.push_back(bFlux);                                               // add in new flux value
                        psi_current[i][j] = front_coeff * (bFlux + dz * q[i] / (2 * mu_n)); // recalculate cell-center flux
                    }
                    bmflux.push_back((1 / (1 - alpha)) * (psi_current[i][j] - alpha * bFlux)); // Append new flux value with alpha!=0
                }
            }
            else if (mu_n < 0)
            {
                fmflux.clear();
                fmflux.push_back(gammas_R[i]);
                for (int j = tot_n_cells - 1; j > -1; j--)
                {
                    region_index = find_region_index(j, n_cells); // determines the spatial region we are in
                    dz = dzs[region_index];
                    Sig_t = Sig_ts[region_index];
                    // Back sweep in z
                    fFlux = fmflux.back();
                    front_coeff = 1 / (1 + Sig_t * dz / (2 * abs(mu_n)));
                    psi_current[i][j] = front_coeff * (fFlux + dz * q[i] / (2 * abs(mu_n)));

                    // This handles when the flux goes negative.
                    // Take the incoming flux to be equal to the previous cell-center flux (alpha=0)
                    if (psi_current[i][j] < 0)
                    {
                        fmflux.erase(std::find(fmflux.begin(), fmflux.end(), fFlux));               // Delete the final flux in fmflux, the value that made us go negative
                        fFlux = psi_current[i][j + 1];                                              // reassign flux for alpha=0 only in this cell
                        fmflux.push_back(fFlux);                                                    // add in new flux value
                        psi_current[i][j] = front_coeff * (fFlux + dz * q[i] / (2 * abs(mu_n))); // recalculate cell-center flux
                    }

                    fmflux.push_back(1 / (1 - alpha) * (psi_current[i][j] - alpha * fFlux)); // Append new flux value with alpha!=0
                }
            }
        }

        // Create vectors from the columns of the angular fluxes to integrate
        for (int m = 0; m < tot_n_cells; m++)
        {
            region_index = find_region_index(m, n_cells); // determines the spatial region we are in
            Sig_s = Sig_ss[region_index];
            for (int k = 0; k < n_mu; k++)
            {
                angular_flux_col.push_back(psi_current[k][m]); // This will loop over a particular column in psi_current to get a vector
            }
            scalar_flux_current.push_back(GQ_integrate(angular_flux_col, w_vec) * 2 * M_PI); // Calculate value of scalar flux at given z-value

            for (int k = 0; k < n_mu; k++) // loop over angles to update scatter source
            {
                q[k] = Sig_s / (4 * M_PI) * scalar_flux_current[m] + S[k]; // Update scatter source --> scalar flux is true over all angles
            }
            angular_flux_col.clear(); // Clear angular flux vector                                            // Clear angular flux vector
        }

        // ----------------------------------------------------------------------
        // Calculate relative error -- use L2 norm
        // ----------------------------------------------------------------------
        if (iteration_num > 0)
        {
            // error = calc_error(scalar_flux_current, scalar_flux_old);

            // reset values
            error_vec.clear();
            running_error = 0;

            // begin L2 norm calculation
            for (int i = 0; i < scalar_flux_current.size(); i++)
            {
                region_index = find_region_index(i, n_cells); // determines the spatial region we are in
                dz = dzs[region_index];
                error_vec.push_back(scalar_flux_old[i] - scalar_flux_current[i]); // take the difference
                running_error += error_vec.back() * error_vec.back() * dz;        // square and multiply by the spacing
            }
            error = sqrt(running_error);
            printf("        Iteration %d: Error = %e\n",
                   iteration_num, error);
        }
        iteration_num += 1;
        if (iteration_num == max_iters)
        {
            printf("MAX ITERATIONS REACHED! "
                   "FINAL ERROR = %e\n",
                   error);
            break;
        }

        // ----------------------------------------------------------------------
        // Clear values and set up for next iteration
        // ----------------------------------------------------------------------
        scalar_flux_old.clear();
        for (int j = 0; j < tot_n_cells; j++)
        {
            scalar_flux_old.push_back(scalar_flux_current[j]);
            for (int k = 0; k < n_mu; k++)
            {
                psi_old[k][j] = psi_current[k][j];
                psi_current[k][j] = 0; // Reset to 0
            }
        }
        scalar_flux_current.clear();
    }

    printf("    Completed processing %d cells.\n", tot_n_cells);

    // ----------------------------------------------------------------------
    // Create output files
    // ----------------------------------------------------------------------

    // Values of the scalar flux will necessarily be reported at the cell CENTER

    printf("    Exporting scalar fluxes to CSV file...\n");

    output_filename = infile_name + "_out.csv";
    std::ofstream outfile(output_filename);

    outfile << "Number of Cells:," << tot_n_cells << "\n";
    outfile << "Final Error:," << error << "\n\n";

    outfile << "z,Flux\n";

    double z_sum = dzs[0] / 2.0;
    for (int j = 0; j < tot_n_cells; j++)
    {
        region_index = find_region_index(j, n_cells); // determines the spatial region we are in
        dz = dzs[region_index];
        outfile << z_sum << "," << scalar_flux_old[j] << "\n";
        z_sum += dz;
    }
    outfile.close();
    printf("    CSV export completed: %s\n\n", output_filename.c_str());

    printf("\nSN2 Code Completed Successfully.\n\n\n");
    return 0;
}
