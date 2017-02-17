//
//  Tests.c
//  quarks EOS
//
//  Created by Clebson Graeff on 2016-06-08.
//  Copyright © 2016 Clebson Graeff. All rights reserved.
//

#include "Tests.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "EOS.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

#include "Constants.h"


void CurvesOfFig2_9Buballa2005(double chem_pot_min_value,
                               double chem_pot_max_value,
                               double mass_guess,
                               double renorm_chem_pot_guess,
                               int    points_number);

void RunTests()
{
#pragma mark Reproduce Fig. 2.9 of M. Buballa, Physics Reports
    // Reproduce Fig. 2.9 of M. Buballa, Physics Reports 407 (2005) 205-376
    // This is a test of the different kinds of transition for mass
    // and density
    if (true)
    {
        printf("\t- Reproduce Fig. 2.9 from  M. Buballa,\n"
               "\t  Physics Reports 407 (2005) 205-376\n");
        SetFilePath("tests/Buballa-Fig1-R/data/");

        double points_number = 1000;
        double chem_pot_min_value = 1.0E-3;
        double chem_pot_max_value = 1000.0;
        double mass_guess = 400.0;
        double renorm_chem_pot_guess = 0.0;

        SetParametersSet("BuballaR_2_GV");

        parameters.model.G_V = 0.0;

        CurvesOfFig2_9Buballa2005 (chem_pot_min_value,
                                   chem_pot_max_value,
                                   mass_guess,
                                   renorm_chem_pot_guess,
                                   points_number);

        parameters.model.G_V = 0.5 * parameters.model.G_S;

        CurvesOfFig2_9Buballa2005 (chem_pot_min_value,
                                   chem_pot_max_value,
                                   mass_guess,
                                   renorm_chem_pot_guess,
                                   points_number);

        parameters.model.G_V = parameters.model.G_S;

        CurvesOfFig2_9Buballa2005 (chem_pot_min_value,
                                   chem_pot_max_value,
                                   mass_guess,
                                   renorm_chem_pot_guess,
                                   points_number);

        // Chiral limit
        parameters.model.bare_mass = 0.0;

        parameters.model.G_V = 0.0;

        CurvesOfFig2_9Buballa2005 (chem_pot_min_value,
                                   chem_pot_max_value,
                                   mass_guess,
                                   renorm_chem_pot_guess,
                                   points_number);

        parameters.model.G_V = 0.5 * parameters.model.G_S;

        CurvesOfFig2_9Buballa2005 (chem_pot_min_value,
                                   chem_pot_max_value,
                                   mass_guess,
                                   renorm_chem_pot_guess,
                                   points_number);

        parameters.model.G_V = parameters.model.G_S;

        CurvesOfFig2_9Buballa2005 (chem_pot_min_value,
                                   chem_pot_max_value,
                                   mass_guess,
                                   renorm_chem_pot_guess,
                                   points_number);

        SetFilePath(NULL);
    }

#pragma mark Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611

    // Reproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611 (1996) 393-408
    // (the figure uses parameters of Set II from the article)
    if (false)
    {
        printf("\tReproduce Fig. 1 (right) from  M. Buballa, Nuclear Physics A 611\n");
        SetFilePath("tests/Buballa-Fig1-R/data/");

        SetParametersSet("Buballa_2");

        double chemical_potential[4] = {0.0, 350.0, 378.5, 410.0};

        double vacuum_mass = VacuumMassDetermination();

        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass,
                                                                       0.0,
                                                                       0.0,
                                                                       0.0);

        for (int i = 0; i < 4; i++){

            double minimum_mass = 0.0;
            double maximum_mass = 1000.0;
            int points_number = 1000;

            double step = (maximum_mass - minimum_mass) / (points_number - 1);

            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);

            double m = 0;

            for (int j = 0; j < points_number; j++) {

                double renorm_chem_pot = chemical_potential[i];

                double fermi_momentum = FermiMomentum(m, renorm_chem_pot);

                double thermodynamic_potential =
                    ThermodynamicPotential(m,
                                           fermi_momentum,
                                           chemical_potential[i],
                                           renorm_chem_pot);

                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output,
                               j,
                               thermodynamic_potential - vacuum_thermodynamic_potential);

                m += step;
            }

            char filename[256];
            sprintf(filename, "Fig1_Buballa_%d.dat", i);

            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential\n",
                               2,
                               mass_vector,
                               output);
        }

        SetFilePath("tests/Buballa-Fig1-R/");
        FILE * log_file = OpenFile("run.log");
        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of vacuum mass, resultin in %f MeV.\n",
                vacuum_mass);
        fprintf(log_file,
                "\tCalculation of the thermodynamic potential as function of mass "
                "(Reproduce Fig. 1 (right) from  M. Buballa, "
                "Nuclear Physics A 611 (1996) 393-408).\n"
                "\tfiles: tests/data/Fig1_Buballa_*.dat\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports

    // Reproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports 407 (2005) 205-376
    // (the figure uses parameters of Set II from the article, with G_V = 0)
    if (false)
    {
        printf("\tReproduce Fig. 2.8 (left) from  M. Buballa, Physics Reports\n");

        SetFilePath("tests/Buballa-Fig2.8-L/data/");

        SetParametersSet("BuballaR_2");

        int points_number = 1000;

        parameters.model.bare_mass = 0; // In the reference, the bare mass was
                                        // zero for this test

        double mass = 0;
        parameters.variables.temperature = 0.0;

        double chemical_potential[4] = {0.0, 300.0, 368.6, 400.0};

        double min_renormalized_chemical_potential = 0.0;
        double max_renormalized_chemical_potential = 1000;

        double min_mass = 0.0;
        double max_mass = 1000.0;

        double renorm_chem_pot_step = Step(min_renormalized_chemical_potential,
                                           max_renormalized_chemical_potential,
                                           points_number);

        double mass_step = Step(min_mass, max_mass, points_number);

        double vacuum_mass = VacuumMassDetermination();

        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass,
                                                                       0.0,
                                                                       0.0,
                                                                       0.0);

        for (int i = 0; i < 4; i++){

            char filename_1[256];
            sprintf(filename_1, "ZeroedRenormalizedChemPotEquation_BR2L_%d.dat", i);

            FILE * f = OpenFile(filename_1);

            renorm_chem_pot_equation_input input;
            input.chemical_potential = chemical_potential[i];
            input.temperature = parameters.variables.temperature;
            input.mass = mass;

            double mu_r = 0;
            while (mu_r < points_number) {
                fprintf(f,
                        "%20.15E\t%20.15E\n",
                        mu_r,
                        ZeroedRenormalizedChemicalPotentialEquation(mu_r,
                                                                    &input));
                mu_r += renorm_chem_pot_step;
            }

            fclose(f);

            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);

            double m = 0;

            for (int j = 0; j < points_number; j++) {

                double renorm_chem_pot = chemical_potential[i];

                double fermi_momentum = FermiMomentum(m, renorm_chem_pot);

                double thermodynamic_potential =
                    ThermodynamicPotential(m,
                                           fermi_momentum,
                                           chemical_potential[i],
                                           renorm_chem_pot);

                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output,
                               j,
                               thermodynamic_potential - vacuum_thermodynamic_potential);

                m += mass_step;
            }

            char filename[256];
            sprintf(filename, "Fig2.8L_BuballaR_%d.dat", i);

            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential\n",
                               2,
                               mass_vector,
                               output);

        }

        SetFilePath("tests/Buballa-Fig2.8-L/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of the hermodynamic potential as function of mass "
                "(Reproduce Fig. 2.8 (left) from  M. Buballa, "
                "Physics Reports 407 (2005) 205-376.\n"
                "\tFiles:tests/data/Fig2.8L_BuballaR_*.dat\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports

    // Reproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports 407 (2005) 205-376
    // (the figure uses parameters of Set II from the article, with G_V = G_S)
    if (false)
    {
        printf("\tReproduce Fig. 2.8 (right) from  M. Buballa, Physics Reports\n");
        SetFilePath("tests/Buballa-Fig2.8-R/data/");

        SetParametersSet("BuballaR_2_GV");

        int points_number = 1000;

        parameters.variables.temperature = 0.0;

        parameters.model.bare_mass = 0; // In the reference, the bare mass was
                                        // zero for this test

        double chemical_potential[4] = {0.0, 430.0, 440.0, 444.3};

        double mass[4] = {1.0, 100, 300, 700};

        double min_renorm_chemical_potential = 0.0;
        double max_renorm_chemical_potential = 1000.0;

        double mu_r_step = Step(min_renorm_chemical_potential,
                                max_renorm_chemical_potential,
                                points_number);

        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++){
                char filename_1[256];
                sprintf(filename_1,
                        "ZeroedRenormalizedChemPotEquation_BR2R_%d_%d.dat",
                        i,
                        j);

                FILE * f = OpenFile(filename_1);

                renorm_chem_pot_equation_input input;
                input.chemical_potential = chemical_potential[i];
                input.mass = mass[j];
                input.temperature = parameters.variables.temperature;

                double mu_r = 0;
                while (mu_r < points_number) {
                    fprintf(f,
                            "%20.15E\t%20.15E\n",
                            mu_r,
                            ZeroedRenormalizedChemicalPotentialEquation(mu_r,
                                                                        &input));
                    mu_r += mu_r_step;
                }

                fclose(f);
            }
        }

        double vacuum_mass = VacuumMassDetermination();

        double vacuum_thermodynamic_potential = ThermodynamicPotential(vacuum_mass,
                                                                       0.0,
                                                                       0.0,
                                                                       0.0);

        for (int i = 0; i < 4; i++){

            double minimum_mass = 0.0;
            double maximum_mass = 1000.0;
            int points_number = 1000;

            double step = (maximum_mass - minimum_mass) / (points_number - 1);

            gsl_vector * mass_vector = gsl_vector_alloc(points_number);
            gsl_vector * output = gsl_vector_alloc(points_number);
            gsl_vector * renorm_chem_pot_vector = gsl_vector_alloc(points_number);

            // Prepare function to be passed to the root finding algorithm
            renorm_chem_pot_equation_input input;
            gsl_function F;
            F.function = &ZeroedRenormalizedChemicalPotentialEquation;
            F.params = &input;

            /*
             * Determine the renormalized chemical potential by root-finding
             * for each value of mass
             */

            double m = 0;
            for (int j = 0; j < points_number; j++) {

                double renorm_chem_pot;

                // Prepare input for ZeroedRenormalizedChemicalPotentialEquation
                input.chemical_potential = chemical_potential[i];
                input.mass = m;

                // If the chemical potential is zero, the solution is zero.
                // Otherwise it needs to be calculated
                if (chemical_potential[i] == 0){
                    renorm_chem_pot = 0;
                }
                else{
                    const double lower_bound = 1.0E-3;
                    const double upper_bound = 1000.0;
                    const double abs_error = 1.0E-5;
                    const double rel_error = 1.0E-5;
                    const double max_iter = 2000;

                    int status = UnidimensionalRootFinder(&F,
                                                          lower_bound,
                                                          upper_bound,
                                                          abs_error,
                                                          rel_error,
                                                          max_iter,
                                                          &renorm_chem_pot);
                    if (status == -1)
                        renorm_chem_pot = 0;
                }

                gsl_vector_set(renorm_chem_pot_vector,
                               j,
                               renorm_chem_pot);

                double thermodynamic_potential =
                    ThermodynamicPotential(m,
                                           chemical_potential[i],
                                           renorm_chem_pot,
                                           parameters.variables.temperature);

                gsl_vector_set(mass_vector, j, m);
                gsl_vector_set(output,
                               j,
                               thermodynamic_potential - vacuum_thermodynamic_potential);

                m += step;
            }

            char filename[256];
            sprintf(filename, "Fig2.8R_BuballaR_%d.dat", i);

            WriteVectorsToFile(filename,
                               "# mass, thermodynamic potential\n",
                               2,
                               mass_vector,
                               output);

            sprintf(filename, "renormalized_chemical_potential_%d.dat", i);
            WriteVectorsToFile(filename,
                               "#mass, renormalized_chemical_potential\n",
                               2,
                               mass_vector,
                               renorm_chem_pot_vector);

            gsl_vector_free(output);
            gsl_vector_free(mass_vector);
            gsl_vector_free(renorm_chem_pot_vector);

        }

        SetFilePath("tests/Buballa-Fig2.8-R/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file,
                "The following tests were executed for %s parameterization:\n",
                parameters.parameters_set_identifier);
        fprintf(log_file,
                "\tCalculation of the thermodynamic potential as function of mass "
                "(Reproduce Fig. 2.8 (right) from  M. Buballa, "
                "Physics Reports 407 (2005) 205-376.\n"
                "\tFiles:tests/data/Fig2.8R_BuballaR_*.dat\n");
        fprintf(log_file,
                "\tZeroed renormalized chemical potential equation "
                "for selected parameters"
                " (as function of renormalized chemical potential)\n");
        fprintf(log_file,
                "\tRenormalized chemical potential (as function of mass)\n");
        fprintf(log_file,
                "\tVacuum mass: %f\n", vacuum_mass);
        fprintf(log_file,
                "\tVacuum thermodynamic potential: %f\n",
                vacuum_thermodynamic_potential);

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }

#pragma mark Reproduce Fig. 2.7 from Buballa Physics Reports
    // TODO: Check this test. We have to make it work for
    // zero density, however we are now running on chemical
    // potential
    if (false)
    {
        printf("\tReproduce Fig. 2.7 from Buballa Physics Reports\n");
        SetFilePath("tests/Buballa-Fig2.7/data/");

        SetParametersSet("BuballaR_2");

        int n_pts = 1000;

        double renorm_chem_pot = 0.0; // would this make barionic_density == 0?
                                      //
        double temperature_min = 0.0;
        double temperature_max = 300.0;
        double temperature_step = Step (temperature_min,
                                        temperature_max,
                                        n_pts);

        gsl_vector * temperature_vector = gsl_vector_alloc(n_pts);
        gsl_vector * mass_vector = gsl_vector_alloc(n_pts);

        double temperature = temperature_min;
        for (int i = 0; i < n_pts; i++){

            parameters.variables.temperature = temperature;

            gap_equation_input params;
            params.renormalized_chemical_potential = 0.0; // maybe this guarantees zero density?
            params.temperature = 0.0;

            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            F.function = &ZeroedGapEquation;
            F.params = &params;

            double mass;
            int status = UnidimensionalRootFinder(&F,
                                                  0.0,
                                                  500.0,
                                                  1.0E-7,
                                                  1.0E-7,
                                                  1000,
                                                  &mass);
            if (status == -1)
                mass = 0;

            gsl_vector_set(temperature_vector, i, temperature);
            gsl_vector_set(mass_vector, i, mass);

            temperature += temperature_step;
        }

        WriteVectorsToFile("mass_temperature.dat",
                           "# Reproduce Fig. 2.7 from Buballa, "
                           "Physics Reports bare_mass not zero\n"
                           "# temperature, mass\n",
                           2,
                           temperature_vector,
                           mass_vector);

        // Repeat for bare_mass = 0
        parameters.model.bare_mass = 0;

        temperature = temperature_min;
        for (int i = 0; i < n_pts; i++){

            parameters.variables.temperature = temperature;

            gap_equation_input input;
            input.renormalized_chemical_potential = renorm_chem_pot;

            // Prepare function to be passed to the root finding algorithm
            gsl_function F;
            F.function = &ZeroedGapEquation;
            F.params = &input;

            // This case will not have solutions beyond T = 220 MeV
            double mass = 0;

            if (temperature < 220.0){
                int status = UnidimensionalRootFinder(&F,
                                                      0.1,
                                                      500.0,
                                                      1.0E-7,
                                                      1.0E-7,
                                                      6000,
                                                      &mass);
                if (status == -1)
                    mass = 0;
            }

            gsl_vector_set(temperature_vector, i, temperature);
            gsl_vector_set(mass_vector, i, mass);

            temperature += temperature_step;
        }

        WriteVectorsToFile("mass_temperature_bare_mass_zero.dat",
                           "# Reproduce Fig. 2.7 from Buballa, "
                           "Physics Reports bare_mass = 0\n"
                           "# temperature, mass\n",
                           2,
                           temperature_vector,
                           mass_vector);

        gsl_vector_free(temperature_vector);
        gsl_vector_free(mass_vector);

        SetFilePath("tests/Buballa-Fig2.7/");
        FILE * log_file = OpenFile("run.log");

        fprintf(log_file, "Reproduce Fig. 2.7 from Buballa Physics Reports\n");

        PrintParametersToFile(log_file);

        fclose(log_file);

        SetFilePath(NULL);
    }
}

void CurvesOfFig2_9Buballa2005(double chem_pot_min_value,
                               double chem_pot_max_value,
                               double mass_guess,
                               double renorm_chem_pot_guess,
                               int points_number)
{
    double chemical_potential = chem_pot_min_value;
    double chemical_potential_step = Step(chem_pot_min_value,
                                          chem_pot_max_value,
                                          points_number);

    gsl_vector * chemical_potential_vector = gsl_vector_alloc(points_number);
    gsl_vector * mass_vector = gsl_vector_alloc(points_number);
    gsl_vector * renorm_chem_pot_vector = gsl_vector_alloc(points_number);
    gsl_vector * barionic_density_vector = gsl_vector_alloc(points_number);

    // Solution of Gap Equation, determination of scalar density
    // the values of these two variables will be updated and reused
    // as guesses by SimultaneousSolution
    double mass = mass_guess;
    double renormalized_chemical_potential = renorm_chem_pot_guess;
    for (int i = 0; i < points_number; i++){

        gsl_vector_set(chemical_potential_vector, i, chemical_potential);

        SimultaneousSolution(parameters.variables.temperature,
                             chemical_potential,
                             mass,
                             renormalized_chemical_potential,
                             &mass,
                             &renormalized_chemical_potential);

        gsl_vector_set(mass_vector, i, mass);
        gsl_vector_set(renorm_chem_pot_vector,
                       i,
                       renormalized_chemical_potential);

        double barionic_density = BarionicDensity(mass,
                                                  renormalized_chemical_potential,
                                                  parameters.variables.temperature);

        gsl_vector_set(barionic_density_vector, i, barionic_density);

        chemical_potential += chemical_potential_step;
    }

    // Use an index to diferentiate all
    // files that will ve written
    static int file_index = 0;

    char filename[256];

    sprintf(filename, "mass_vs_chem_pot_%d", file_index);
    WriteVectorsToFile (filename,
                        "# chemical potential (MeV), mass (MeV) \n",
                        2,
                        chemical_potential_vector,
                        mass_vector);

    sprintf(filename, "bar_dens_vs_chem_pot_%d.dat", file_index);
    WriteVectorsToFile (filename,
                        "# chemical potential (MeV), barionic density (fm^{-3})\n",
                        2,
                        chemical_potential_vector,
                        barionic_density_vector);

    file_index++;

    gsl_vector_free(chemical_potential_vector);
    gsl_vector_free(mass_vector);
    gsl_vector_free(renorm_chem_pot_vector);
    gsl_vector_free(barionic_density_vector);

    return;
}
