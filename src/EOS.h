//
//  EOS.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef FiniteTemperatureEOS_h
#define FiniteTemperatureEOS_h

// Parameters for the determination of
// entropy by integration
typedef struct _EntropyParameters{
    double lower_limit;
    double upper_limit;
    double abs_error;
    double rel_error;
    int max_sub_interval;
    int integration_key;
} EntropyParameters;

typedef struct _SimultaneousSolutionParameters{
    int max_iter;
    double mass_guess; // (MeV)
    double renorm_chem_pot_guess; //(MeV)
    double abs_error;
    double rel_error;

    UnidimensionalRootFindingParameters zero_mass_case;
    EntropyParameters entropy;

} SimultaneousSolutionParameters;

void SimultaneousSolution(double temperature,
                          double chemical_potential,
                          double mass_guess,
						  double renor_chem_pot_guess,
                          double *return_mass,
                          double *return_renorm_chem_pot);

double BarionicDensity(double mass,
                       double renormalized_chemical_potential,
                       double temperature);

double ScalarDensity(double temperature,
                     double mass,
                     double renormalized_chemical_potential);

double FermiMomentum(double mass,
                     double renormalized_chemical_potential);

double ThermodynamicPotential(double mass,
                              double chemical_potential,
                              double renormalized_chemical_potential,
                              double temperature);

double ThermodynamicPotentialFreeGasContribution(double mass,
                                                 double chemical_potential,
                                                 double renormalized_chemical_potential,
                                                 double temperature);

double Pressure(double regularized_thermodynamic_potential,
                double temperature);

double EnergyDensity(double regularized_thermodynamic_potential,
                     double chemical_potential,
                     double barionic_density,
                     double temperature,
                     double entropy);

double Entropy(double mass, double temperature, double renormalized_chemical_potential);
double EntropyIntegrand(double momentum, void * parameters);

typedef struct _gap_equation_input{
    double renormalized_chemical_potential;
	double temperature;
} gap_equation_input;

typedef struct _renorm_chem_pot_equation_input{
    double chemical_potential;
    double mass;
	double temperature;
} renorm_chem_pot_equation_input;

typedef struct _multi_dim_root_params{
    double chemical_potential;
    double temperature;
} multi_dim_root_params;

double ZeroedGapEquation(double mass,
                         void * params);

double ZeroedRenormalizedChemicalPotentialEquation(double  renorm_chemical_potential,
                                                   void   *params);

double VacuumMassDetermination();
double VacuumMassEquation(double mass, void * input);

#endif /* FiniteTemperatureEOS_h */
