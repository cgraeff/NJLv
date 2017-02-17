//
//  Parameters.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Constants.h"
#include "Parameters.h"
#include "CommandlineOptions.h"

// Chosen parameter set (globally accessible)
Parameters parameters;
static Parameters parameters_sets_list[256] = {0};
static int parameters_sets_list_count = 0;

Parameters NewCopyOfParametersSetFromTemplate();
void AppendParametersSetToList(Parameters a_set);

void ParametersSetup(void)
{
    // See header file for units and other relevant
    // information about the parameters below.

    // START DECLARATION OF PARAMETERS SETS:
    //
    //  For declaration of parameters sets:
    //      - copy from the default one,
    //      - overwrite values as desired,
    //      - append to the list of parameterizations
    //
    // The first declared parameters set will be treated as the standard
    // case. To run others, use the "-p IDENTIFIER" commandline option, without
    // quotation marks) where IDENTIFIER is the content of the variable
    // parameters_set_identifier in the definitions bellow.

    Parameters p;

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "Buballa_1";

    p.model.parameters_set_origin = "Set 1 from M. Buballa, "
                                    "Nucl. Phys. A 611 (1996) 393-408";

    p.model.bare_mass = 0.0;    // MeV
    p.model.cutoff = 650.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.14 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.0;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "Buballa_2";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, "
                                    "Nucl. Phys. A 611 (1996) 393-408";

    p.model.bare_mass = 0.0;    // MeV
    p.model.cutoff = 600.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.45 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.0;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "Buballa_3";

    p.model.parameters_set_origin = "Set 3 from M. Buballa, "
                                    "Nucl. Phys. A 611 (1996) 393-408";

    p.model.bare_mass = 0.0;    // MeV
    p.model.cutoff = 570.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.84 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.0;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics "
                                    "Reports 407 (2005) 205-376 with G_V = 0";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.0;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics "
                                    "Reports 407 (2005) 205-376 with G_V = G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV_0.25";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports"
                                    "407 (2005) 205-376 with G_V = 0.25 * G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.25 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV_0.35";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.35 * G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.35 * p.model.G_S;

    AppendParametersSetToList(p);


    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV_0.45";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.45 * G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.45 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV_0.55";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.55 * G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.55 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV_0.65";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.65 * G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.65 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "BuballaR_2_GV_0.75";

    p.model.parameters_set_origin = "Set 2 from M. Buballa, Physics Reports "
                                    "407 (2005) 205-376 with G_V = 0.75 * G_S";

    p.model.bare_mass = 5.6;    // MeV
    p.model.cutoff = 587.9;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.75 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "PCP-0.0";

    p.model.parameters_set_origin = "Set from PRD 94 094001, 2016";

    p.model.bare_mass = 5.1;    // MeV
    p.model.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.11 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.0 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "PCP-0.1";

    p.model.parameters_set_origin = "Set from PRD 94 094001, 2016";

    p.model.bare_mass = 5.1;    // MeV
    p.model.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.11 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.1 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "PCP-0.2";

    p.model.parameters_set_origin = "Set from PRD 94 094001, 2016";

    p.model.bare_mass = 5.1;    // MeV
    p.model.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.11 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.2 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "PCP-0.3";

    p.model.parameters_set_origin = "Set from PRD 94 094001, 2016";

    p.model.bare_mass = 5.1;    // MeV
    p.model.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.11 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.3 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "PCP-0.4";

    p.model.parameters_set_origin = "Set from PRD 94 094001, 2016";

    p.model.bare_mass = 5.1;    // MeV
    p.model.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.11 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.4 * p.model.G_S;

    AppendParametersSetToList(p);

    ////

    p = NewCopyOfParametersSetFromTemplate();
    p.parameters_set_identifier = "PCP-0.5";

    p.model.parameters_set_origin = "Set from PRD 94 094001, 2016";

    p.model.bare_mass = 5.1;    // MeV
    p.model.cutoff = 648.0;     // MeV

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0)
    // corrects the dimension
    p.model.G_S = 2.11 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.5 * p.model.G_S;

    AppendParametersSetToList(p);

    // END DECLARATION OF PARAMETERS SETS

    // Verify that the set identifiers are unique
    for (int i = 0; i < parameters_sets_list_count; i++){
        for (int j = 0; j < parameters_sets_list_count; j++){

            if (i == j)
                continue;

            if(!strcasecmp(parameters_sets_list[i].parameters_set_identifier,
                           parameters_sets_list[j].parameters_set_identifier)){
                printf("Two parameters sets share the \"%s\" identifier,"
                       " but it should be unique.\n",
                       parameters_sets_list[i].parameters_set_identifier);
                exit(EXIT_FAILURE);
            }
        }
    }

    // If asked to, print parameters sets and exit
    if (options.list_available_parameterizations){
        for (int i = 0; i < parameters_sets_list_count; i++){
            printf("Parameters set %s\n"
                   "\tOrigin: %s\n",
                   parameters_sets_list[i].parameters_set_identifier,
                   parameters_sets_list[i].model.parameters_set_origin);
        }
        exit(EXIT_SUCCESS);
    }

    // If there isn't at least one set, exit
    if (parameters_sets_list_count == 0){
        printf("There are not parameters set declared. Declare at least one set.\n");
        exit(EXIT_FAILURE);
    }
}

void SetParametersSet(char * parameters_set_identifier)
{
    // If the identifier is null, return default case
    if (parameters_set_identifier == NULL){
        parameters = parameters_sets_list[0];
        return;
    }

    for (int i = 0; i < parameters_sets_list_count; i++){

        Parameters p = parameters_sets_list[i];

        if (!strcasecmp(p.parameters_set_identifier, parameters_set_identifier)){
            parameters = p;
            return;
        }
    }

    printf("Parameters set %s unrecognized.\n"
           "Use -l to list available parameterizations.\n",
           parameters_set_identifier);
    exit(EXIT_FAILURE);
}

Parameters NewCopyOfParametersSetFromTemplate()
{
    Parameters p;

    p.parameters_set_identifier = "Template";

    p.model.parameters_set_origin = "A standard parameters set. With theory parameters "
                                    "from Set 2 of M. Buballa, Physics Reports 407"
                                    "(2005) 205-376";
    p.model.bare_mass = 5.6;
    p.model.cutoff = 587.9;

    // The code expects [G_S] = fm^2, the pow(CONST_HBAR_C, 2.0) corrects the dimension
    p.model.G_S = 2.44 * pow(CONST_HBAR_C / p.model.cutoff, 2.0);
    p.model.G_V = 0.0;

    p.variables.points_number = 1000;
    p.variables.min_value = 1.0E-3;
    p.variables.max_value = 2000.0;
    p.variables.temperature = 0.0; // (MeV)

    // Low lower_bound but not zero, as it may be problematic if bare_mass == 0
    // upper_bound near the value of the nucleon mass.

    p.vacuum_mass_determination.lower_bound = 1.0E-3;
    p.vacuum_mass_determination.upper_bound = 1000.0;
    p.vacuum_mass_determination.abs_error = 1.0E-5;
    p.vacuum_mass_determination.rel_error = 1.0E-5;
    p.vacuum_mass_determination.max_iterations = 2000;

    // In the following, the guesses are used for the first iteration.
    // Any subsequent iteration uses as guesses the values of mass
    // and renormalized chemical potential from the previous iteration.
    // This shall work well for stepping of a control variable.
    p.simultaneous_solution.max_iter = 4000;
    p.simultaneous_solution.mass_guess = 300.0; // (MeV)
    p.simultaneous_solution.renorm_chem_pot_guess = 400.0; //(MeV)
    p.simultaneous_solution.abs_error = 1.0E-5;
    p.simultaneous_solution.rel_error = 1.0E-5;
    p.simultaneous_solution.zero_mass_case.renorm_chem_pot_lower_bound = 0.0;
    p.simultaneous_solution.zero_mass_case.renorm_chem_pot_upper_bound = 3000;

    p.fermi_dirac_integrals.lower_limit = 0.0;
    p.fermi_dirac_integrals.upper_limit = p.model.cutoff;
    p.fermi_dirac_integrals.max_interval_num = 8000;
    p.fermi_dirac_integrals.integration_key = GSL_INTEG_GAUSS61;
    p.fermi_dirac_integrals.max_sub_interval = 8000;
    p.fermi_dirac_integrals.abs_error = 1.0E-10;
    p.fermi_dirac_integrals.rel_error = 1.0E-10;

    p.therm_pot_free_gas_integral.lower_limit = 0.0;
    p.therm_pot_free_gas_integral.upper_limit = p.model.cutoff;
    p.therm_pot_free_gas_integral.max_interval_num = 8000;
    p.therm_pot_free_gas_integral.integration_key = GSL_INTEG_GAUSS61;
    p.therm_pot_free_gas_integral.max_sub_interval = 8000;
    p.therm_pot_free_gas_integral.abs_error = 1.0E-10;
    p.therm_pot_free_gas_integral.rel_error = 1.0E-10;

    p.simultaneous_solution.entropy.lower_limit = 0.0;
    p.simultaneous_solution.entropy.upper_limit = parameters.model.cutoff;
    p.simultaneous_solution.entropy.abs_error = 1.0E-3;
    p.simultaneous_solution.entropy.rel_error = 1.0E-3;
    p.simultaneous_solution.entropy.max_sub_interval = 1000;
    p.simultaneous_solution.entropy.integration_key = GSL_INTEG_GAUSS61;

    return p;
}

void AppendParametersSetToList(Parameters a_set)
{
    // Append to list of parameterizations
    parameters_sets_list[parameters_sets_list_count] = a_set;
    parameters_sets_list_count++;
}

void PrintParametersToFile(FILE * file)
{
    fprintf(file, "\n--- PARAMETERS ---\n");
    fprintf(file, "set_identifier = %s\n\n", parameters.parameters_set_identifier);

    fprintf(file, "set_origin = %s\n", parameters.model.parameters_set_origin);
    fprintf(file, "G_S = %f\n", parameters.model.G_S);
    fprintf(file, "G_V = %f\n", parameters.model.G_V);
    fprintf(file, "cutoff = %f\n", parameters.model.cutoff);
    fprintf(file, "bare_mass = %f\n", parameters.model.bare_mass);

    fprintf(file, "points_number = %d\n", parameters.variables.points_number);
    fprintf(file, "minimum_density = %f\n", parameters.variables.min_value);
    fprintf(file, "maximum_density = %f\n\n", parameters.variables.max_value);

    fprintf(file,
            "vac_mass_det_max_iterations = %d\n",
            parameters.vacuum_mass_determination.max_iterations);
    fprintf(file,
            "vac_mass_det_lower_bound = %f\n",
            parameters.vacuum_mass_determination.lower_bound);
    fprintf(file,
            "vac_mass_det_upper_bound = %f\n",
            parameters.vacuum_mass_determination.upper_bound);
    fprintf(file,
            "vac_mass_det_abs_error = %f\n",
            parameters.vacuum_mass_determination.abs_error);
    fprintf(file,
            "vac_mass_det_rel_error = %f\n\n",
            parameters.vacuum_mass_determination.rel_error);

    fprintf(file,
            "temperature = %f\n\n",
            parameters.variables.temperature);
    fprintf(file,
            "simultaneous_solution.max_iter = %d\n",
            parameters.simultaneous_solution.max_iter);
    fprintf(file,
            "simultaneous_solution.mass_guess = %f\n",
            parameters.simultaneous_solution.mass_guess);
    fprintf(file,
            "simultaneous_solution.renorm_chem_pot_guess = %f\n",
            parameters.simultaneous_solution.renorm_chem_pot_guess);
    fprintf(file,
            "simultaneous_solution.abs_error = %f\n",
            parameters.simultaneous_solution.abs_error);
    fprintf(file,
            "simultaneous_solution.rel_error = %f\n\n",
            parameters.simultaneous_solution.rel_error);

    fprintf(file,
            "fermi_dirac_integrals.lower_limit = %f\n",
            parameters.fermi_dirac_integrals.lower_limit);
    fprintf(file,
            "fermi_dirac_integrals.upper_limit = %f\n",
            parameters.fermi_dirac_integrals.upper_limit);
    fprintf(file,
            "fermi_dirac_integrals.max_interval_num = %d\n",
            parameters.fermi_dirac_integrals.max_interval_num);
    fprintf(file,
            "fermi_dirac_integrals.integration_key = %d\n",
            parameters.fermi_dirac_integrals.integration_key);
    fprintf(file,
            "fermi_dirac_integrals.max_sub_interval = %d\n",
            parameters.fermi_dirac_integrals.max_sub_interval);
    fprintf(file,
            "fermi_dirac_integrals.abs_error = %f\n",
            parameters.fermi_dirac_integrals.abs_error);
    fprintf(file,
            "fermi_dirac_integrals.rel_error = %f\n\n",
            parameters.fermi_dirac_integrals.rel_error);

    return;
}

