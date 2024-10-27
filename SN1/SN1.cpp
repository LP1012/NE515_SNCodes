// To run code:
// g++ -I./src SN1.cpp ./src/*.cpp -o SN1 && ./SN1
// Need to find a way to make VSCode do this...

// This may (or may not) be faster:
// g++ -I./src SN1.cpp ./src/*.cpp -o SN1 -DAE_CPU=AE_INTEL -mavx2 -mfma -DAE_OS=AE_POSIX -pthread && ./SN1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <cmath>

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

double calc_error(const std::vector<double> &current_flux, const std::vector<double> &previous_flux)
// Function determines the current error by calculating the relative error between scalar flux vectors using the 2-norm
{
    double cfsum = 0, pfsum = 0;
    for (int k = 0; k < current_flux.size(); k++)
    {
        cfsum += current_flux[k] * current_flux[k];
        pfsum += previous_flux[k] * previous_flux[k];
    }
    double cfnorm = sqrt(cfsum);
    double pfnorm = sqrt(pfsum);

    double error = std::abs(cfnorm - pfnorm) / pfnorm;

    return error;
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

int main(int argc, char const *argv[])
{
    printf("\n***************************************************\n");
    printf("                 SN1 CODE EXECUTION\n");
    printf("***************************************************\n\n");
    double Sig_s = 0.5, L = 5, Sig_t = 1.0, Sig_f = 0.0; // given material data
    double tol = 1e-6;                                   // convergence tolerance
    double dz;                                           // allocate memory for later
    int n_cells[4] = {10, 20, 50, 100};                  // Define number of cells in space
    double alpha = 0.5;                                  // specifies diamond differencing
    int max_iters = 100;                                 // Safety measure
    alglib::ae_int_t n_mu = 8;                           // Number of quadrature nodes in mu
    std::vector<double> fmflux;
    std::vector<double> bmflux;
    string output_filename;

    printf("\nRun Specifications:\n"
           "    Sigma_s: %.4f            Sigma_f: %.4f\n"
           "    Sigma_t: %.4f            Length (L): %.1f\n"
           "    Alpha: %.4f              Max Iterations: %d\n"
           "    Error Tolerance: %e      Quadrature Nodes: %ld\n\n",
           Sig_s, Sig_f, Sig_t, L, alpha, max_iters, tol, n_mu);

    // Variables to store the output
    alglib::ae_int_t info;          // To hold the status code
    alglib::real_1d_array mus_init; // To store nodes
    alglib::real_1d_array w;        // To store weights
    std::vector<double> w_vec;      // Store weights as vector
    double mus[n_mu] = {};          // Store discrete ordinates in array

    // Call the Gauss-Legendre quadrature generator
    printf("Calculating GQ weights and discrete ordinates... ");
    alglib::gqgenerategausslegendre(n_mu, info, mus_init, w);
    printf("Done.\n"
           "Casting GQ data to a vector (weights) and array (scattering angle)...");

    // Cast weights into a vector and mu's into an array
    for (int i = 0; i < n_mu; i++)
    {
        w_vec.push_back(w[i]);
        mus[i] = mus_init[i];
    }
    printf("Done.\n"
           "Initializing variables and iterating over cell sizes...\n");
    std::cout << "\nGQ Weights:  w = " << w.tostring(4);
    std::cout << "\nGQ Nodes:   mu = " << mus_init.tostring(4) << "\n";

    double mu_n, z_i, front_coeff;        // Allocate memory for variables in loops
    double fFlux, bFlux;                  // Corresponding to "front flux" and "back flux", respectively
    std::vector<double> angular_flux_col; // Vector to store column of angular flux values (corresponding to the same point in space)

    // Define BCs
    double gamma_L = 5000000;
    double gamma_R = 0;

    // ----------------------------------------------------------------------
    // Begin looping over number of cells in spatial coordinates
    // ----------------------------------------------------------------------

    printf("\n---------------------------------------------\n"
           "Begin looping over cell numbers...\n"
           "---------------------------------------------\n\n");
    for (int I : n_cells)
    {
        printf("Number of cells = %d\n", I);
        printf("    Set spatial step size... ");
        dz = L / (double)I;
        printf("dz = %.4f\n", dz);

        // Create 2D arrays to hold the scalar flux at different mu and z values
        // mu is varied in rows, while z is varied by column
        printf("    Initializing angular flux arrays...");
        double psi_current[n_mu][I] = {};
        double psi_old[n_mu][I] = {};
        printf("Done.\n");

        // Must initialize q because it is used in the first iteration.
        // When compared to the notes, this will be q^(n-1), and code will simply make a direct update
        // instead of creating a new vector.
        printf("    Intializing scattering source terms... ");
        double q[I] = {};
        double S[I] = {}; // Initialize a source vector as well

        for (int i = 0; i < I; i++)
        {
            q[i] = 0; // Initialize with all 0's
            S[i] = 0; // CHANGE IF NEEDED FOR FUTURE APPLICATIONS
        }
        printf("Done.\n");

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
                    bFlux = gamma_L;
                    for (int j = 0; j < I; j++)
                    {
                        // Forward sweep in z
                        front_coeff = 1 / (1 + Sig_t * dz / (2 * mu_n));
                        psi_current[i][j] = front_coeff * (bFlux + dz * q[j] / (2 * mu_n));
                        bFlux = (1 / (1 - alpha)) * (psi_current[i][j] - alpha * bFlux);
                    }
                }
                else if (mu_n < 0)
                {
                    fFlux = gamma_R;
                    // printf("Beginning backward sweep...\n");
                    for (int j = I - 1; j > -1; j--)
                    {
                        // Back sweep in z
                        front_coeff = 1 / (1 + Sig_t * dz / (2 * abs(mu_n)));
                        // printf("DEBUG: front_coeff = %f\n",front_coeff);
                        psi_current[i][j] = front_coeff * (fFlux + dz * q[j] / (2 * abs(mu_n)));
                        fFlux = 1 / (1 - alpha) * (psi_current[i][j] / alpha * fFlux);
                    }
                    // printf("Backward sweep completed.\n\n");
                }
            }

            // Create vectors from the columns of the angular fluxes to integrate
            for (int m = 0; m < I; m++)
            {
                for (int k = 0; k < n_mu; k++)
                {
                    angular_flux_col.push_back(psi_current[k][m]); // This will loop over a particular column in psi_current to get a vector
                }
                scalar_flux_current.push_back(GQ_integrate(angular_flux_col, w_vec));   // Calculate value of scalar flux at given z-value
                q[m] = Sig_s / (4 * M_PI) * scalar_flux_current[m] + S[m] / (4 / M_PI); // Update scatter source
                angular_flux_col.clear();                                               // Clear angular flux vector
            }

            // ----------------------------------------------------------------------
            // Calculate relative error
            // ----------------------------------------------------------------------
            if (iteration_num > 0)
            {
                error = calc_error(scalar_flux_current, scalar_flux_old);
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
            for (int j = 0; j < I; j++)
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

        printf("    Completed processing %d cells.\n", I);

        // ----------------------------------------------------------------------
        // Create output files
        // ----------------------------------------------------------------------
        
        // Values of the scalar flux will necessarily be reported at the cell CENTER

        printf("    Exporting scalar fluxes to CSV file...\n");

        output_filename = to_string(I) + "_out.csv";
        std::ofstream outfile(output_filename);

        outfile << "Left Boundary:," << gamma_L << "\n";
        outfile << "Right Boundary:," << gamma_R << "\n";
        outfile << "Number of Cells:," << I << "\n";
        outfile << "Final Error:," << error << "\n\n";

        outfile << "z,Flux\n";

        for (int j = 0; j < I; j++)
        {
            outfile << dz * (double)(j + 1) - dz / 2.0 << "," << scalar_flux_old[j] << "\n";
        }
        outfile.close();
        printf("    CSV export completed: %s\n\n", output_filename.c_str());
    }

    printf("\nSN1 Code Completed Successfully.\n\n\n");
    return 0;
}
