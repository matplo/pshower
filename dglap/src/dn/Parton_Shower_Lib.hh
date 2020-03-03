#ifndef Parton_Shower_Lib
#define Parton_Shower_Lib
#include <gsl/gsl_rng.h>

#define as 0.1187
#define MZ 91.87
#define PI 3.14159265359
#define DIPOLE_INITIAL_CAPACITY 200000

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

typedef struct {
  long int size;      // slots used so far
  long int capacity;  // total available slots
  long int *emissions_Label;
  long int *Flavor;
  long int *Dead;
  double *Polar;
  double *Azimuth;
  double *Efrac;
  double *ParentSplitAngle;
} PSEmissionsList;

int PSemissions_double_capacity_if_full(PSEmissionsList *emissions);

void PSemissions_init(PSEmissionsList *emissions);
void PSemissions_append(PSEmissionsList *emissions, double *value, long int PartonFlavor, long int PartonDead, long int Label);
void PSemissions_edit(PSEmissionsList *emissions, long int index, double *value, long int PartonFlavor, long int PartonDead, long int Label);
void RetrieveParton(PSEmissionsList *emissions, long int TheChosen, double *Momentum, long int *Flavor, long int *Death, long int *Label);
void PSemissions_free(PSEmissionsList *emissions);
int  WriteNewPSEmission(PSEmissionsList *emissions, double *Momentum, long int PartonFlavor,
				   long int PartonDead, long int *EmissionsLabel);

////////////////////////////////////////////////////////////////


//We solve for RF in the equation t(RI,RF) = - Int[as(m)/m,{m,Q*RI,Q*RF}]/PI, RF < RI
double ComputeSplitAngle(double Q, double RI, double t, double  NF, double CA);

//We solve for KT in the equation t(RI,RF) = - Int[as(m)/m,{m,Q,KT}]/PI, RF < RI
double ComputeKTRunningCoupling(double Q, double RI, double t, double  NF, double CA);

double Alphas(double MU, double NF, double TheCA);

double AlphasIntegral(double Q, double RI, double RF, double NF, double TheCA);

double LambertW(const double z);

////////////////////////////////////////////////////////////////

double ComputeVirtuality(double Efrac, double Q, double Z, double SplitAngle);
double ComputeKT(double Efrac, double Q, double Z, double SplitAngle);
double MinZFunction(double Efrac, double Q, double PhysSplitAngle, double MinQ, double p, double q );

////////////////////////////////////////////////////////////////

/*
We note the vales in params:
*(params) = Minimum Energy Fraction or Zcut;
*(params+1) = CA;
*(params+2) = NF;
*(params+3) = CF;
*(params+4) = Jet Energy;
*(params+5) = Jet Initial Opening Angle;
*(params+6) = Minimum Mass Scale;
*/

void FindAZwithVeto(long int CurrentFlavor, double Efrac, double *Daughter1, double *Daughter2,
			   long int *AssignFlavor, double *Z, double *params,
			   const gsl_rng *r);

void SplitSelectedPartonWithVeto(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
					long int *TheChosen, double *params, const gsl_rng *r, long int *EndEvent, long int Mode);

void SplitSelectedParton(PSEmissionsList *emissions,  PSEmissionsList *DaughterEmissions, double *Angle, double *t,
				long int *TheChosen, double *params, const gsl_rng *r, long int *EndEvent, long int Mode,
				double *InvertPg, double *InvertPq, long int NumInverts);

void BuildSplitFunctionInversion(double *InvertPg, double *InvertPq, long int NumInverts, double *params, const gsl_rng *r);

void DGLAPDownToAngle(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t,
			     double *params, double RF, double *ParentMomentum, long int *ParentFlavor, long int *ParentLabel, long int *TheChosen,
			     double *InvertPg, double *InvertPq, long int NumInverts, long int *EndEvent, long int Mode, const gsl_rng *r);

void DGLAPByVetoDownToAngle(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *Angle, double *t,
				   double *params, double RF, double *ParentMomentum, long int *ParentFlavor, long int *ParentLabel,
				   long int *TheChosen, long int *EndEvent, long int Mode, const gsl_rng *r);

void WriteToDiskHistgramsFRAG(char filename[1000], long int TotalNumEvents, long int NumofRadii,
				     long int NumBins, double **BINS, double *Radii, double TheMinZBook,
				     long int Flavor, double *params, long int LogBin);

void ReInitialize(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, long int Flavor, double *t);

void LoadShowerParams(char filename[1000], double *params, long int *NumRadii, double *JetRadii);

void UpdateEmissionList(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis,
			       double *ParentMomentum, long int *ParentFlavor, long int *ParentLabel, long int TheChosen);

////////////////////////////////////////////////////////////////

#endif
