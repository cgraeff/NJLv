//
//  main.c
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_vector.h>

#include "libdatafun/libdatafun.h"

#include "CommandlineOptions.h"
#include "Parameters.h"
#include "Constants.h"
#include "Tests.h"
#include "EOS.h"

#include "Loop.h"

int main(int argc, char * argv[])
{
    CommandlineOptionsParse(argc, argv);
    ParametersSetup();

    if (options.tests){
        RunTests();
        exit(EXIT_SUCCESS);
    }

    // If option -p is used, set parameters set accordingly,
    // otherwise, use default set
    SetParametersSet(options.parameterization);

    // If the temperature was chosen using
    // commandline options, use it
    // (-1.0 is a place holder value)
    if (options.temp != -1.0)
        parameters.variables.temperature = options.temp;

    if (parameters.variables.temperature >= 0){
        SolveFiniteTemperatureEOS();
    }
    else{
        printf("Values of temperature must be non-negative.\n");
        printf("(%f was provided).\n",
               parameters.variables.temperature);
        exit(EXIT_FAILURE);
    }

    return 0;
}

