#include "jet_fragmentation.hh"
#include <math.h>

void test(double **LeadingJetSpectra)
{
  ;
}

double **doubleParray(long int n, long int k)
{
  double **p = (double **) malloc( (n + 1) * sizeof(double *));
  for(long int i=0; i < n + 1; i++)
  {
    p[i] = (double *) malloc( ( k + 1 ) * sizeof(double ));
    for (long int j=0; j < k + 1; j++)
      p[i][j] = 0;
  }
  return p;
}

void JetFragmentation(PSEmissionsList *emissions, PSEmissionsList *DaughterEmissions, double *CurrentWTAaxis, double *JetRadii,
		      long int NumRadii, double *params, long int NumEbins,
		      double **LeadingJetSpectra, double **EventWideSpectra, double **EventWideSpectraLogBin,
		      long int NumEvent, double *t, double *InvertPg, double *InvertPq, long int NumInverts,
		      long int Mode, const gsl_rng *r){

  long int i,j,k,Ebin,LogEbin,TheChosen,CurrentFlavor,CurrentLabel,Death,EndEvent, CurrentMaxEbin;
  double Angle,RI,RF,FindLogEbinStart, MinZ =*(params);
  double ParentMomentum[4],Momentum[4];
  long int ParentFlavor,ParentLabel;
  /*
    Momentum[4]={Polar,Azimuth,Efrac,SplitAngle}
   */

  *(t) = 0;
  i = 0;
  j = 0;

  //Before we start evolution, we log the histograms of the zero-th
  //order emission
  if((emissions->size)!=1){printf("What, more than one emission to start?\n");}

  while(i<NumRadii){
    RF = *(JetRadii+i);

    DGLAPDownToAngle(emissions, DaughterEmissions, CurrentWTAaxis, &Angle, t,
		     params, RF, ParentMomentum, &ParentFlavor, &ParentLabel, &TheChosen,
		     InvertPg, InvertPq, NumInverts, &EndEvent, Mode, r);

    //Now we update the event wide histograms
    CurrentMaxEbin = 0;
    for(j=0; j < (emissions->size); j++){
      RetrieveParton(emissions, j, Momentum, &CurrentFlavor, &Death, &CurrentLabel);

      Ebin = (int)floor(NumEbins * Momentum[2] );
      if(Ebin>CurrentMaxEbin){
	CurrentMaxEbin = Ebin;
      };
      FindLogEbinStart = NumEbins*(1 - log(Momentum[2])/log(MinZBook) );
      if(FindLogEbinStart > 0){
	LogEbin = (int)floor(FindLogEbinStart);
      }else{
	LogEbin = 0;
      }

      EventWideSpectra[i][Ebin] = EventWideSpectra[i][Ebin] + 1;
      EventWideSpectraLogBin[i][LogEbin] = EventWideSpectraLogBin[i][LogEbin] + 1;
    };//for

    LeadingJetSpectra[i][CurrentMaxEbin] = LeadingJetSpectra[i][CurrentMaxEbin] + 1;

    //The very last splitting was not folded into the emissions.
    //we do that now, and resume evolution.
    UpdateEmissionList(emissions, DaughterEmissions, CurrentWTAaxis,
		       ParentMomentum, &ParentFlavor, &ParentLabel, TheChosen );

    i=i+1;
  };
  return;
}
