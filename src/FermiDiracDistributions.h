//
//  FermiDiracDistributions.h
//  NJLv
//
//  Created by Clebson Graeff on 2017-02-14.
//  Copyright © 2017 Clebson Graeff. All rights reserved.
//

#ifndef FermiDiracDistributions_h
#define FermiDiracDistributions_h

double FermiDiracDistributionForParticles(double energy,
                                          double chemical_potential,
                                          double temperature);

double FermiDiracDistributionForAntiparticles(double energy,
                                              double chemical_potential,
                                              double temperature);

double FermiDiracDistributionFromDensityIntegral(double temperature,
                                                 double mass,
                                                 double renormalized_chemical_potential);

double FermiDiracDistributionIntegralFromScalarDensity(double temperature,
                                                       double mass,
                                                       double renorm_chemical_potential);

#endif /* FermiDiracDistributions_h */
