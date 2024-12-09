// To run code:
// g++ -I./src SN1.cpp ./src/*.cpp -o SN1 && ./SN1

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

double calc_error(const std::vector<double> current_flux, const std::vector<double> previous_flux, double dz)
// Function determines the current error by calculating the relative error between scalar flux vectors using the discrete L2 norm
// Numerical integration is implemented using the midpoint rule
{
    double cfsum = 0, pfsum = 0;
    for (int k = 0; k < current_flux.size(); k++)
    {
        cfsum += current_flux[k] * current_flux[k];
        pfsum += previous_flux[k] * previous_flux[k];
    }
    double cfnorm = sqrt(dz * cfsum);
    double pfnorm = sqrt(dz * pfsum);

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
        return 0;
    }
}

double Q_mms(double sigt, double mu, double z)
// function returns the value of the source term evaluated at a particular
// location (z) and angle (mu) for the manufacture solution
{
    double mus = 5.0 * pow(mu, 3.0) - 3.0 * mu;
    double z1 = 0.8 - 8.0 / 25.0 * z;
    double z2 = -4.0 / 25.0 * pow(z, 2.0) + 0.8 * z;

    double Qm = 0.5 * mus * (mu * z1 + sigt * z2);
    return Qm;
}

// ----------------------------------------------------------------------

int main(int argc, char const *argv[])
{
    printf("\n***************************************************\n");
    printf("                 SN1 CODE EXECUTION\n");
    printf("***************************************************\n\n");
    double Sig_s = 1.0, L = 5.0, Sig_t = 1.0, Sig_f = 0.0; // given material data
    double tol = 1e-6;                                     // convergence tolerance
    double dz;                                             // allocate memory for later
    int n_cells[2] = {20, 100};                            // Define number of cells in space
    double alpha = 0.5;                                    // specifies diamond differencing
    int max_iters = 1000;                                  // Safety measure
    alglib::ae_int_t n_mu = 8;                             // Number of quadrature nodes in mu
    std::vector<double> fmflux;                            // allocate memory for storing cell fluxes -- "forward marching flux"
    std::vector<double> bmflux;                            // allocate memory for storing cell fluxes -- "backward marching flux"
    std::vector<double> RH_boundary_fluxes;                // store fluxes for RH albedo condition
    double albedo = 0.5;                                   // albedo value
    string output_filename;

    bool zero_dirichlet_bool = false;

    printf("\nRun Specifications:\n"
           "    Sigma_s: %.4f               Sigma_f: %.4f\n"
           "    Sigma_t: %.4f               Length (L): %.1f\n"
           "    Alpha: %.4f                 Max Iterations: %d\n"
           "    Quadrature Nodes: %ld           Error Tolerance: %.1e\n\n",
           Sig_s, Sig_f, Sig_t, L, alpha, max_iters, n_mu, tol);

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
        printf("    Homogeneous Dirichlet conditions specified. Scaling factor calculation ignored.\n"
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
        double j_left = 0;
        for (int k = 0; k < n_mu; k++)
        {
            gammas_L[k] *= (1 / BC_scaling_factor_left);

            j_left += gammas_L[k] * mus[k] * w_vec[k];

            // *** UNCOMMENT THIS FOR NONZERO RH BOUNDARY CONDITION ***
            // gammas_R[k] *= (1/BC_scaling_factor_right);
            RH_boundary_fluxes.push_back(0.0); // assigns RH albedo condtion with zeros for first iteration
        }
        cout << "j_left = " << 2.0 * M_PI * j_left << endl;
        printf("Done.\n");
    }

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
        printf("    Initializing angular flux arrays... ");
        double psi_current[n_mu][I] = {};
        double psi_old[n_mu][I] = {};
        printf("Done.\n");

        // Must initialize q because it is used in the first iteration.
        // When compared to the notes, this will be q^(n-1), and code will simply make a direct update
        // instead of creating a new vector.
        printf("    Intializing scattering source terms... ");
        // Because of angular dependence, both must be 2D arrays.
        // Same structure as the angular flux (psi) values
        double q[n_mu][I] = {};
        double S[n_mu][I] = {};

        double z_current, z_next; // variables to be used in quadrature
        for (int i = 0; i < n_mu; i++)
        // loop over angles
        {
            mu_n = mus[i];
            for (int j = 0; j < I; j++)
            // loop over cells
            {
                // implement a trapezoid quadrature rule for cell-averaged source
                z_current = (double)j * dz;
                z_next = double(j + 1) * dz;

                // The following implements a quadrature rule for a source term.
                // Since this code was not intended to solve for a source other than zero,
                // the structure is left as a placeholder.
                S[i][j] = 1 / 2.0 * (0 + 0); // delta_z cancels for volume-average
                q[i][j] = S[i][j];           // for the first iteration, we assume a 0 scalar flux everywhere
            }
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
                    bmflux.clear();
                    bmflux.push_back(gammas_L[i]);
                    for (int j = 0; j < I; j++)
                    {
                        // Forward sweep in z
                        bFlux = bmflux.back();
                        front_coeff = 1 / (1 + Sig_t * dz / (2 * mu_n));
                        psi_current[i][j] = front_coeff * (bFlux + dz * q[i][j] / (2 * mu_n));

                        // This handles when the flux goes negative.
                        // Take the incoming flux to be equal to the previous cell-center flux (alpha=0)
                        if (psi_current[i][j] < 0)
                        {
                            bmflux.erase(std::find(bmflux.begin(), bmflux.end(), bFlux));          // Delete the final flux in bmflux, the value that made us go negative
                            bFlux = psi_current[i][j - 1];                                         // reassign flux for alpha=0 only in this cell
                            bmflux.push_back(bFlux);                                               // add in new flux value
                            psi_current[i][j] = front_coeff * (bFlux + dz * q[i][j] / (2 * mu_n)); // recalculate cell-center flux
                        }
                        bmflux.push_back((1 / (1 - alpha)) * (psi_current[i][j] - alpha * bFlux)); // Append new flux value with alpha!=0

                        if (j == 0)
                        {
                            RH_boundary_fluxes[i] = bmflux.back();
                        }
                    }
                }
                else if (mu_n < 0)
                {
                    fmflux.clear();

                    RH_boundary_fluxes[i] = albedo * RH_boundary_fluxes[n_mu - (i + 1)];
                    fmflux.push_back(albedo * RH_boundary_fluxes[n_mu - (i + 1)]);

                    cout << "DEBUG: mu_n = " << mu_n << endl;
                    cout << "DEBUG: Albedo flux for i = " << i << ": " << fmflux.back() << endl;

                    // fmflux.push_back(gammas_R[i]);
                    for (int j = I - 1; j > -1; j--)
                    {
                        // Back sweep in z
                        fFlux = fmflux.back();
                        front_coeff = 1 / (1 + Sig_t * dz / (2 * abs(mu_n)));
                        psi_current[i][j] = front_coeff * (fFlux + dz * q[i][j] / (2 * abs(mu_n)));

                        // This handles when the flux goes negative.
                        // Take the incoming flux to be equal to the previous cell-center flux (alpha=0)
                        if (psi_current[i][j] < 0)
                        {
                            fmflux.erase(std::find(fmflux.begin(), fmflux.end(), fFlux));               // Delete the final flux in fmflux, the value that made us go negative
                            fFlux = psi_current[i][j + 1];                                              // reassign flux for alpha=0 only in this cell
                            fmflux.push_back(fFlux);                                                    // add in new flux value
                            psi_current[i][j] = front_coeff * (fFlux + dz * q[i][j] / (2 * abs(mu_n))); // recalculate cell-center flux
                        }

                        fmflux.push_back(1 / (1 - alpha) * (psi_current[i][j] - alpha * fFlux)); // Append new flux value with alpha!=0
                    }
                }
            }

            // Create vectors from the columns of the angular fluxes to integrate
            for (int m = 0; m < I; m++)
            {
                for (int k = 0; k < n_mu; k++)
                {
                    angular_flux_col.push_back(psi_current[k][m]); // This will loop over a particular column in psi_current to get a vector
                }
                scalar_flux_current.push_back(GQ_integrate(angular_flux_col, w_vec) * 2 * M_PI); // Calculate value of scalar flux at given z-value

                for (int k = 0; k < n_mu; k++) // loop over angles to update scatter source
                {
                    q[k][m] = Sig_s / (4 * M_PI) * scalar_flux_current[m] + S[k][m]; // Update scatter source --> scalar flux is true over all angles
                }
                angular_flux_col.clear(); // Clear angular flux vector                                            // Clear angular flux vector
            }

            // ----------------------------------------------------------------------
            // Calculate relative error
            // ----------------------------------------------------------------------
            if (iteration_num > 0)
            {
                error = calc_error(scalar_flux_current, scalar_flux_old, dz);
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

        output_filename = to_string(I) + "_case5_albedo-0.5_out.csv";
        std::ofstream outfile(output_filename);

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
