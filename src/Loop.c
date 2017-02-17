//
//  Loop.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <gsl/gsl_vector.h>

#include "libdatafun/headers/AuxiliaryFunctions.h"
#include "libdatafun/headers/Derivatives.h"

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "EOS.h"

int SolveFiniteTemperatureEOS(){

    // Print name of parametrization
    if (options.verbose){
        printf("Calculation performed with %s parameters set\n"
               "\ttemperature: %4.2f (MeV)\n"
               "\tmin chemical potential: %4.2E (MeV)\n"
               "\tmax chemical potential: %4.2E (MeV)\n",
               parameters.parameters_set_identifier,
               parameters.variables.temperature,
               parameters.variables.min_value,
               parameters.variables.max_value);
    }

    // Vectors to store results
    int points_number = parameters.variables.points_number;

    gsl_vector * barionic_density_vector = gsl_vector_alloc(points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(points_number);
    gsl_vector * chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * renormalized_chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * thermodynamic_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * pressure_vector = gsl_vector_alloc(points_number);
    gsl_vector * entropy_vector = gsl_vector_alloc(points_number);
    gsl_vector * energy_density_vector = gsl_vector_alloc(points_number);


    // Vacuum mass determination
    if (options.verbose)
        printf("Determining the vacuum mass and bag constant ...\n");

    double vacuum_mass = VacuumMassDetermination();

    double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass,
                                                                   0.0,
                                                                   0.0,
                                                                   0.0);

    double bag_constant =
        ThermodynamicPotential(parameters.model.bare_mass, 0.0, 0.0, 0.0)
        - vacuum_thermodynamic_potential;

    if (options.verbose){
        printf("\tVacuum mass: %f\n", vacuum_mass);
        printf("\tBag constant: %f\n", bag_constant);
    }

    // Main loop (on barionic density)
    VariableParameters loop = parameters.variables;

    double chemical_potential = loop.min_value;
    double chemical_potential_step = Step(loop.min_value,
                                          loop.max_value,
                                          loop.points_number);

    if (options.verbose)
        printf("Solving gap equation and equations of state ...\n");

    // Get guesses from parameters
    double mass = parameters.simultaneous_solution.mass_guess;
    double renormalized_chemical_potential =
        parameters.simultaneous_solution.renorm_chem_pot_guess;

    for (int i = 0; i < points_number; i++){

        if (options.verbose){

            double progress = (chemical_potential - loop.min_value)
                              / (loop.max_value - loop.min_value);

            printf("\r\tProgress: %4.1f%%", progress * 100.0);
            fflush(stdout);
        }

        gsl_vector_set(chemical_potential_vector, i, chemical_potential);

        // Solution of Gap Equation, determination of scalar density
        // The function will update the variables for mass and
        // renormalized_chemical_potential and the new values will be used
        // as guesses for the next iteration
        SimultaneousSolution(parameters.variables.temperature,
                             chemical_potential,
                             mass,
                             renormalized_chemical_potential,
                             &mass,
                             &renormalized_chemical_potential);

        gsl_vector_set(mass_vector, i, mass);
        gsl_vector_set(renormalized_chemical_potential_vector,
                       i,
                       renormalized_chemical_potential);

        double barionic_density = BarionicDensity(mass,
                                                  renormalized_chemical_potential,
                                                  parameters.variables.temperature);

        gsl_vector_set(barionic_density_vector, i, barionic_density);

        double thermodynamic_potential =
            ThermodynamicPotential(mass,
                                   chemical_potential,
                                   renormalized_chemical_potential,
                                   parameters.variables.temperature);

        double regularized_thermodynamic_potential = thermodynamic_potential
                                                     - vacuum_thermodynamic_potential;

        // Just the regularized thermodynamic potential will be saved
        // as it's the only one that is used in other calculations
        gsl_vector_set(thermodynamic_potential_vector,
                       i,
                       regularized_thermodynamic_potential);

        // Determination of pressure
        double pressure = Pressure(regularized_thermodynamic_potential,
                                   parameters.variables.temperature);

        gsl_vector_set(pressure_vector, i, pressure);

        // If temperature == 0, we skip the entropy calculation
        // leaving it undefined
        double entropy = 0;
        if (parameters.variables.temperature != 0){
            entropy = Entropy(mass,
                              parameters.variables.temperature,
                              renormalized_chemical_potential);

            gsl_vector_set(entropy_vector, i, entropy);
        }

        // Determination of energy density
        double energy_density = EnergyDensity(regularized_thermodynamic_potential,
                                              chemical_potential,
                                              barionic_density,
                                              parameters.variables.temperature,
                                              entropy);

        gsl_vector_set(energy_density_vector, i, energy_density);

        chemical_potential += chemical_potential_step;
    }
    if (options.verbose)
        printf("\n"); // As print inside the loop doesn't use new line, we need one now

    // Calculate energy per particle
    gsl_vector * energy_density_per_particle_vector =
        VectorNewVectorFromDivisionElementByElement(energy_density_vector,
                                                    barionic_density_vector);

    // Save results
    if (options.verbose)
        printf("Saving results ...\n");

    if (options.dirs)
        SetFilePath("output/IR/data/");

    WriteVectorsToFile ("barionic_density.dat",
                        "# chemical potential (MeV), barionic density (fm^{-3}) \n",
                        2,
                        chemical_potential_vector,
                        barionic_density_vector);

    WriteVectorsToFile("mass.dat",
                       "# chemical potential (MeV), mass (MeV)\n",
                       2,
                       chemical_potential_vector,
                       mass_vector);

    WriteVectorsToFile ("renormalized_chemical_potential.dat",
                        "# chemical potential (MeV), "
                        "renormalized chemical potential (MeV)\n",
                        2,
                        chemical_potential_vector,
                        renormalized_chemical_potential_vector);

    WriteVectorsToFile ("thermodynamic_potential.dat",
                        "# chemical potential (MeV), "
                        "thermodynamic potential (MeV/fm^3) \n",
                        2,
                        chemical_potential_vector,
                        thermodynamic_potential_vector);

    WriteVectorsToFile ("pressure.dat",
                        "# chemical potential (MeV), pressure (MeV/fm^3) \n",
                        2,
                        chemical_potential_vector,
                        pressure_vector);

    WriteVectorsToFile ("energy_density.dat",
                        "# chemical potential (MeV), energy density (MeV/fm^3) \n",
                        2,
                        chemical_potential_vector,
                        energy_density_vector);

    if (options.dirs)
        SetFilePath("output/EOS/data/");

    WriteVectorsToFile("pressure.dat",
                       "# barionic density, pressure \n",
                       2,
                       barionic_density_vector,
                       pressure_vector);

    WriteVectorsToFile("energy_density.dat",
                       "# barionic density, energy density \n",
                       2,
                       barionic_density_vector,
                       energy_density_vector);

    WriteVectorsToFile("energy_density_per_particle.dat",
                       "# barionic density, energy density per particle \n",
                       2,
                       barionic_density_vector,
                       energy_density_per_particle_vector);

    SetFilePath(NULL);

    // Clean up
    gsl_vector_free(barionic_density_vector);
    gsl_vector_free(mass_vector);
    gsl_vector_free(thermodynamic_potential_vector);
    gsl_vector_free(pressure_vector);
    gsl_vector_free(energy_density_vector);
    gsl_vector_free(energy_density_per_particle_vector);


    if (options.verbose)
        printf("Done!\n");

    return 0;
}
