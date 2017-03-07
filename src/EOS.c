//
//  EOS.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>

#include "libdatafun/libdatafun.h"

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"
#include "EOS.h"
#include "FermiDiracDistributions.h"

// Types for internal data passing
typedef struct _multi_dim_gap_eq_param {
    double chemical_potential;
    double temperature;
} multi_dim_gap_eq_param;

typedef struct _therm_pot_free_gas_contrib_params{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} therm_pot_free_gas_contrib_params;

typedef struct _entropy_integrand_parameters{
    double mass;
    double temperature;
    double renormalized_chemical_potential;
} entropy_integrand_parameters;

// Internal helper functions
int MultiDimensionalRootFinderHelperFunction(const gsl_vector *x,
                                             void             *p,
                                             gsl_vector       *return_values);

double ZeroMassSpecialCaseHelperFunction(double x,
                                         void  *par);

double ThermodynamicPotentialFreeGasContributionIntegrand(double momentum,
                                                          void * parameters);

double F0(double mass, double momentum);
double F_E(double mass, double momentum);

/*  Solution for gap equation, beta equilibrium and charge neutrality:
 *  Prototype:
 *          void SimultaneousSolution(double  temperature,
 *                                    double  chemical_potential
 *                                    double *return_mass,
 *                                    double *return_renormalized_chemical_potential);
 *
 *  Purpose:
 *      This function will take parameters like initial guess, errors, tolerance
 *      and call the rootfinding functions in a proper way. It exists to handle
 *      mappings (see below) and to handle the special case of mass = 0 (also
 *      below). It will return the mass which solves the gap function and the
 *      renormalized chemical potential.
 *
 *  Mappings:
 *      The variables for the root finding are assumed to cover the range
 *      (-\infty, +\infty), but that is not the case for the variables
 *      that we are trying to solve. Here both the mass $m$ and the renormalized
 *      chemical potential are such that
 *          $m \in [0, +\infty)$,
 *          $\mu_R \in [0, +\infty)$.
 *      To solve that, we use the mappings:
 *          $m = x^2$,
 *          $\mu_R = y^2$.
 *      The initial guesses must be transformed by inverting the relations
 *      above.
 *
 *  Zero mass special case:
 *      The zero mass case is important as most calculations will be
 *      performed at this particular case, which due to characteristics
 *      of the multidimensional root-finding algorithm, may be problematic to
 *      solve (it works most of the time, but sometimes calculations result in NaNs).
 *      This is due to problems in the calculation of derivatives of
 *      the function with respect to mass which arise from low variability
 *      of the function near zero mass.
 *
 *      This case, however is not the one reached at the start of calculations.
 *      In addition to that, once it is reached, all subsequent calculations are
 *      performed at approximatelly zero mass.
 *
 *      We take these characteristics into account and do the zero mass case
 *      with a special path, where we just assume mass == 0 (this effectivelly
 *      reduces the dimension of the system). This will avoid
 *      any calculation of potentially problematic derivatives. The special
 *      path is triggered by the condition
 *          mass < mass_tolerance
 *      which must be true. The tolerance should be adjusted in Constants.h.
 */

void SimultaneousSolution(double  temperature,
                          double  chemical_potential,
                          double  mass_guess,
                          double  renorm_chem_pot_guess,
                          double *return_mass,
                          double *return_renorm_chem_pot)
{
    // Set up parameters to be passed to helper function
    multi_dim_root_params p;
    p.chemical_potential = chemical_potential;
    p.temperature = temperature;

    // Check for zero mass special case. As mass != 0 is the
    // case that appears first, it is implemented first.
    if (mass_guess > ZERO_MASS_TOL){

        // Set dimension (number of equations|variables to solve|find)
        const int dimension = 2;

        gsl_multiroot_function f;
        f.f = &MultiDimensionalRootFinderHelperFunction;
        f.n = dimension;
        f.params = (void *)&p;

        gsl_vector * initial_guess = gsl_vector_alloc(dimension);
        gsl_vector * return_results = gsl_vector_alloc(dimension);

        gsl_vector_set(initial_guess, 0, sqrt(mass_guess));

        gsl_vector_set(initial_guess, 1, sqrt(renorm_chem_pot_guess));

        int status = MultidimensionalRootFinder(dimension,
                                                &f,
                                                initial_guess,
                                                parameters.simultaneous_solution.abs_error,
                                                parameters.simultaneous_solution.rel_error,
                                                parameters.simultaneous_solution.max_iter,
                                                return_results);

        if (status != 0){
            printf("Something is wrong with the rootfinding.\n");
            abort();
        }

        // Save results in return variables,
        // taking care of the mappinps
        *return_mass = pow(gsl_vector_get(return_results, 0), 2.0);
        *return_renorm_chem_pot = pow(gsl_vector_get(return_results, 1), 2.0);

        // Free vectors
        gsl_vector_free(initial_guess);
        gsl_vector_free(return_results);

        return;
    }
    else{ // Handle special case: Zero mass case

        gsl_function F;
        F.function = &ZeroMassSpecialCaseHelperFunction;
        F.params = &p;

        // Set root bounds observing the mappings
        UnidimensionalRootFindingParameters p =
            parameters.simultaneous_solution.zero_mass_case;

        p.lower_bound = sqrt(p.lower_bound);
        p.upper_bound = sqrt(p.upper_bound);

        // As we are left with just one variable and one equation to solve,
        // now an one-dimensional algorithm may be employed. Otherwise,
        // the dimension ought to be decreased by one an the multidimensional
        // employed again.
        double return_result;

        int status = UnidimensionalRootFinder(&F,
                                              p,
                                              &return_result);
        if (status != 0){
            printf("\nBounds do not straddle root.\n");
            abort();
        }

        // Save results in return variables,
        // taking care of the mappings
        *return_mass = 0.0;
        *return_renorm_chem_pot = pow(return_result, 2.0);

        return;
    }
}

double ZeroMassSpecialCaseHelperFunction(double  x,
                                         void   *par)
{
    const int dimension = 2;

    gsl_vector * input_values = gsl_vector_alloc(dimension);
    gsl_vector * return_values = gsl_vector_alloc(dimension);

    // Set mass = 0, which is our special case
    gsl_vector_set(input_values, 0, 0);

    // Pass value selected by the root finding routine
    gsl_vector_set(input_values, 1, x);

    MultiDimensionalRootFinderHelperFunction(input_values, par, return_values);

    double return_value = gsl_vector_get(return_values, 1);

    gsl_vector_free(input_values);
    gsl_vector_free(return_values);

    return return_value;
}

int MultiDimensionalRootFinderHelperFunction(const gsl_vector   *x,
                                             void               *params,
                                             gsl_vector         *return_values)
{
    multi_dim_root_params *p = (multi_dim_root_params *)params;

    // The parameters will be passed to the subfunctions Zeroed*

    // Mappings:
    //      The variables for the root finding are assumed to cover the range
    //      (-\infty, +\infty), but that is not the case for the variables
    //      that we are trying to solve.
    //      To solve that, we use the mappings:
    //          $m = x^2$
    //          $\mu_R = y^2$
    //      The initial guesses must be transformed by inverting the relations
    //      above

    const double mass = pow(gsl_vector_get(x, 0), 2.0);
    const double renormalized_chemical_potential = pow(gsl_vector_get(x,1), 2.0);

    gap_equation_input gap_input;
    gap_input.renormalized_chemical_potential = renormalized_chemical_potential;
    gap_input.temperature = p->temperature;

    double zeroed_gap_eq = ZeroedGapEquation(mass, &gap_input);

    renorm_chem_pot_equation_input mu_r_input;
    mu_r_input.chemical_potential = p->chemical_potential;
    mu_r_input.mass = mass;
    mu_r_input.temperature = p->temperature;

    double zeroed_chem_pot_eq =
        ZeroedRenormalizedChemicalPotentialEquation(renormalized_chemical_potential,
                                                    &mu_r_input);

    gsl_vector_set(return_values, 0, zeroed_gap_eq);
    gsl_vector_set(return_values, 1, zeroed_chem_pot_eq);

    return GSL_SUCCESS;
}

double ZeroedGapEquation(double mass,
                         void * params)
{
    gap_equation_input * p = (gap_equation_input *)params;

    double term = 2.0 * parameters.model.G_S * CONST_HBAR_C
                  * ScalarDensity(p->temperature,
                                  mass,
                                  p->renormalized_chemical_potential);

    return mass - parameters.model.bare_mass + term;
}

double ZeroedRenormalizedChemicalPotentialEquation(double renormalized_chemical_potential,
                                                   void * params)
{
    renorm_chem_pot_equation_input * p = (renorm_chem_pot_equation_input *)params;

    double term = 2.0 * parameters.model.G_V * NUM_COLORS * CONST_HBAR_C
                  * BarionicDensity(p->mass,
                                    renormalized_chemical_potential,
                                    p->temperature);

    return renormalized_chemical_potential - p->chemical_potential + term;
}

double BarionicDensity(double mass,
                       double renormalized_chemical_potential,
                       double temperature)
{
    double constant = NUM_FLAVORS / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum_3rd_power =
            pow(FermiMomentum(mass, renormalized_chemical_potential), 3.0);

        return constant * fermi_momentum_3rd_power / 3.0;
    }

    double integral =
        FermiDiracDistributionFromDensityIntegral(temperature,
                                                  mass,
                                                  renormalized_chemical_potential);

    return constant * integral;
}

double ThermodynamicPotential(double mass,
                              double chemical_potential,
                              double renormalized_chemical_potential,
                              double temperature)
{
    double first_term =
        ThermodynamicPotentialFreeGasContribution(mass,
                                                  chemical_potential,
                                                  renormalized_chemical_potential,
                                                  temperature);

    double second_term = pow(mass - parameters.model.bare_mass, 2.0)
                         / (4.0 * parameters.model.G_S * CONST_HBAR_C);

    // If G_V == 0, we have to avoid a division by zero
    double third_term = 0.0;
    if (parameters.model.G_V != 0)
        third_term = pow(chemical_potential - renormalized_chemical_potential, 2.0)
        / (4.0 * parameters.model.G_V * CONST_HBAR_C);

    return first_term + second_term - third_term;
}

double ThermodynamicPotentialFreeGasContribution(double mass,
                                                 double chemical_potential,
                                                 double renormalized_chemical_potential,
                                                 double temperature)
{
    double constant = - NUM_FLAVORS * NUM_COLORS * pow(CONST_HBAR_C, -3.0)
                        / pow(M_PI, 2.0);

    if (temperature == 0){

        double fermi_momentum = FermiMomentum(mass, renormalized_chemical_potential);

        double F_diff = F_E(mass, parameters.model.cutoff) - F_E(mass, fermi_momentum);

        return constant
               * (F_diff
                  + renormalized_chemical_potential * pow(fermi_momentum, 3.0) / 3.0);
    }

    therm_pot_free_gas_contrib_params p;
    p.mass = mass;
    p.renormalized_chemical_potential = renormalized_chemical_potential;
    p.temperature = temperature;

    gsl_function F;
    F.function = &ThermodynamicPotentialFreeGasContributionIntegrand;
    F.params = &p;

    double integral = OnedimensionalIntegrator(&F,
                                               parameters.therm_pot_free_gas_integral);

    return constant * integral;
}

double ThermodynamicPotentialFreeGasContributionIntegrand(double momentum,
                                                          void * parameters)
{
    therm_pot_free_gas_contrib_params * p =
        (therm_pot_free_gas_contrib_params *)parameters;

    double T = p->temperature;
    double mu_r = p->renormalized_chemical_potential;

    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    // From docs:
    // If x is nearly zero, then the common expression log(1 + x) will
    // not be able to produce accurate results, as most (or all) of the
    // information in x will be lost by addition.  Instead, use
    // log1p(x) to perform the same computation without undue
    // loss of accuracy.
    double first_term = T * log1p(exp(-(energy - mu_r)/T));
    double second_term = T * log1p(exp(-(energy + mu_r)/T));

    return pow(momentum, 2.0) * (energy + first_term + second_term);
}

double Entropy(double mass, double temperature, double renormalized_chemical_potential)
{
    EntropyParameters params = parameters.simultaneous_solution.entropy;

    entropy_integrand_parameters p;
    p.mass = mass;
    p.temperature = temperature;
    p.renormalized_chemical_potential = renormalized_chemical_potential;

    gsl_function F;
    F.function = &EntropyIntegrand;
    F.params = &p;

    gsl_integration_workspace * workspace =
        gsl_integration_workspace_alloc(params.max_sub_interval);

    double integral;
    double abserr;
    gsl_integration_qag(&F,
                        params.lower_limit,
                        params.upper_limit,
                        params.abs_error,
                        params.rel_error,
                        params.max_sub_interval,
                        params.integration_key,
                        workspace,
                        &integral,
                        &abserr);

    gsl_integration_workspace_free(workspace);

    return NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -3.0) * integral
           / pow(M_PI, 2.0);
}

double entropy_integrand_expression(double temperature,
                                    double energy,
                                    double renormalized_chemical_potential)
{
    double arg = - (energy - renormalized_chemical_potential)
                   / temperature;

    return log1p(exp(arg)) - arg * exp(arg) / (1.0 + exp(arg));
}

double EntropyIntegrand(double momentum, void * parameters)
{
    entropy_integrand_parameters * p = (entropy_integrand_parameters *)parameters;

    double energy = sqrt(pow(momentum, 2.0) + pow(p->mass, 2.0));

    double first_term =
        entropy_integrand_expression(p->temperature,
                                     energy,
                                     p->renormalized_chemical_potential);

    double second_term =
        entropy_integrand_expression(p->temperature,
                                     energy,
                                     -p->renormalized_chemical_potential);

    return pow(momentum, 2.0) * (first_term + second_term);
}

double Pressure(double regularized_thermodynamic_potential, double temperature)
{
    return - regularized_thermodynamic_potential;
}

double EnergyDensity(double regularized_thermodynamic_potential,
                     double chemical_potential,
                     double barionic_density,
                     double temperature,
                     double entropy)
{
    if (temperature == 0){
            return regularized_thermodynamic_potential
                   + NUM_COLORS * chemical_potential * barionic_density;
    }

    return regularized_thermodynamic_potential
           + temperature * entropy
           + chemical_potential * NUM_COLORS * barionic_density;
}

double VacuumMassDetermination()
{
    // Prepare function to be passed to the root finding algorithm.
    // No parameters are needed.
    gsl_function F;
    F.function = &VacuumMassEquation;

    double root;
    int status =
        UnidimensionalRootFinder(&F,
                                 parameters.vacuum_mass_determination,
                                 &root);

    if (status == -1){
        return 0;
    }

    return root;
}

double VacuumMassEquation(double mass, void * input)
{
    double F_diff = F0(mass, parameters.model.cutoff) - F0(mass, 0.0);
    double term = 2.0 * NUM_COLORS * NUM_FLAVORS * pow(CONST_HBAR_C, -2.0)
                  * parameters.model.G_S * mass * F_diff
                  / pow(M_PI, 2.0);

    return mass - parameters.model.bare_mass - term;
}

double ScalarDensity(double temperature,
                     double mass,
                     double renormalized_chemical_potential)
{
    if (mass == 0){
        return 0;
    }

    double constant = - NUM_COLORS * NUM_FLAVORS
                        / (pow(M_PI, 2.0) * pow(CONST_HBAR_C, 3.0));

    if (temperature == 0){

        double fermi_momentum = FermiMomentum (mass, renormalized_chemical_potential);

        return constant * mass * (F0(mass, parameters.model.cutoff)
                                  - F0(mass, fermi_momentum));
    }

    double integral =
        FermiDiracDistributionIntegralFromScalarDensity(temperature,
                                                        mass,
                                                        renormalized_chemical_potential);

    return constant * mass * integral;
}

double F0(double mass, double momentum)
{
    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

    return (1.0 / 2.0) * (momentum * E - pow(mass, 2.0) * log((momentum + E) / mass));
}

double F_E(double mass, double momentum)
{
    if (mass == 0)
        return pow(momentum, 4.0) / 4.0;

    double E = sqrt(pow(mass, 2.0) + pow(momentum, 2.0));

    return (momentum * pow(E, 3.0)
            - 0.5 * pow(mass, 2.0) * momentum * E
            - 0.5 * pow(mass, 4.0) * log ((momentum + E) / mass))
           / 4.0;
}

double FermiMomentum(double mass, double renormalized_chemical_potential)
{
    if (pow(renormalized_chemical_potential, 2.0) < pow(mass, 2.0))
        return 0.0;

    return sqrt(pow(renormalized_chemical_potential, 2.0) - pow(mass, 2.0));
}
