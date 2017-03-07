//
//  Parameters.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef Parameters_h
#define Parameters_h

#include <stdbool.h>
#include <stdio.h>

#include "libdatafun/libdatafun.h"

#include "EOS.h"
#include "FermiDiracDistributions.h"
#include "Loop.h"

// Parameterization variables
typedef struct _ModelParameters{
    char * parameters_set_origin;   // Where the set was taken from
    double G_S;                     // scalar-isoscalar coupling (fm^2)
    double G_V;                     // vector-isoscalar coupling (fm^2)
    double cutoff;                  // \Lambda (MeV)
    double bare_mass;               // (MeV)
} ModelParameters;

// Parameters for the variable for which we are
// performing the calculation of the EOS
typedef struct _VariableParameters{
    int points_number; // Number of points in which the above range will be divide into
    double min_value;  // (fm^-3)
    double max_value;  // (fm^-3)
    double temperature; // (MeV)
} VariableParameters;

typedef struct _parameters
{
    char * parameters_set_identifier;

    VariableParameters variables;
    ModelParameters model;

    UnidimensionalRootFindingParameters vacuum_mass_determination;
    SimultaneousSolutionParameters simultaneous_solution;

	IntegratorParameters fermi_dirac_integrals;
	IntegratorParameters therm_pot_free_gas_integral;

} Parameters;

extern Parameters parameters;

void ParametersSetup(void);
void SetParametersSet(char * parameters_set_identifier);

void PrintParametersToFile(FILE * file);

#endif /* Parameters_h */
