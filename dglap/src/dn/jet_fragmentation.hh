#ifndef DGLAP_DN_JET_FRAGMENTATION_HH
#define DGLAP_DN_JET_FRAGMENTATION_HH
#include "Parton_Shower_Lib.hh"

#define MinZBook 0.00001

void test(double **LeadingJetSpectra);

double **doubleParray(long int n);

void JetFragmentation(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *JetRadii,
		      long int NumRadii, double *params, long int NumEbins,
		      double **LeadingJetSpectra, double **EventWideSpectra, double **EventWideSpectraLogBin,
		      long int NumEvent, double *t, double *InvertPg, double *InvertPq, long int NumInverts,
		      long int Mode, const gsl_rng *r);

#endif
