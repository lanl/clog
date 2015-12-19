/* 

0D burn code.  (Very simple). 

Used for testing ideas, magnitude of effects, etc. 

Written by B. Albright June, 2008.



*/   


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>   // For math prototypes
#include <stdint.h> // For fixed width integer types
#include <limits.h> // For integer limits
#include <float.h>  // For floating point limits
#include <stdarg.h> // For va_list, va_start, va_end

#define BEGIN_PRIMITIVE do
#define END_PRIMITIVE   while(0)

// The following macros give a provide a simple logging capability. Due
// to the way they work, usage needs double parenthesis. That is:
//
//   ERROR(("Could not allocate %i bytes", req));
//
// will print the following message to the log:
//
//   Error at src/module/file.c(34):
//           Could not allocate 45 bytes
//
// Note: Error messages are abortive but MESSAGE and WARNING are not

#define CHECKPOINT() BEGIN_PRIMITIVE {                     \
 print_log( "%s(%i): Checkpoint\n", __FILE__, __LINE__ ); \
} END_PRIMITIVE

#define MESSAGE(args) BEGIN_PRIMITIVE {         \
 print_log( "%s(%i): ", __FILE__, __LINE__ );  \
 print_log args;                               \
 print_log( "\n" );                            \
} END_PRIMITIVE

#define WARNING(args) BEGIN_PRIMITIVE {                      \
 print_log( "Warning at %s(%i):\n\t", __FILE__, __LINE__ ); \
 print_log args;                                            \
 print_log( "\n" );                                         \
} END_PRIMITIVE

#define ERROR(args) BEGIN_PRIMITIVE {                      \
 print_log( "Error at %s(%i):\n\t", __FILE__, __LINE__ ); \
 print_log args;                                          \
 print_log( "\n" );                                       \
 exit(1);                                                 \
} END_PRIMITIVE

#define ALLOCATE(A,LEN,TYPE)                                              \
 if ( !((A)=(TYPE *)malloc((size_t)(LEN)*sizeof(TYPE))) )                \
   ERROR(("Cannot allocate.")); 


// Index within a 2D array of values with NI total value

#define INDEX( I, J, NI ) ( (J)*(NI) + (I) )



#define SEARCH_TABLE(X,Y,XVALS,NX) {                                      \
 int il=0, ih=NX, im;                                                    \
 do {                                                                    \
   im=(ih+il)>>1;                                                        \
   if ( X>XVALS[im] ) il=im; else ih=im;                                 \
 } while ( ih-il>1 );                                                    \
 Y=il;                                                                   \
};                        

// Number of elements in ENDF DT reaction rate tables
#define NUMDT 99

// Linear interpolation of reaction rate table
#define COMPUTE_REACTION_RATE(T,REACTION_RATE)                            \
 BEGIN_PRIMITIVE {                                                       \
   int jl;                                                               \
   double tl, th, invdt;                                                 \
   SEARCH_TABLE(T,jl,tval,NUMDT);                                        \
   tl=tval[jl]; th=tval[jl+1]; invdt=1./(th-tl);                         \
   REACTION_RATE=invdt*((T-tl)*svbar[jl+1]+(th-T)*svbar[jl]);            \
 } END_PRIMITIVE; 

// Numerical constants
#define TWO_THIRDS (0.666666667)

// Physical constants
#define RAD_CONST  (0.0137e16)           /* ergs keV^-4           */ 
#define KBOLTZ     (1.602e-9)            /* ergs/keV              */ 
#define CVAC       (2.9979e10)           /* cm/s                  */ 
#define NA         (6.02e23)             /* dimensionless         */ 
#define EALPHA     (3.51e3*KBOLTZ)       /* ergs - alpha energy   */ 
#define MASS_AMU   (1.660538783e-24)     /* atomic mass unit in g */
#define MASS_E     (9.10938215e-28)      /* electron mass in g    */
#define MASS_D     (2.013553*MASS_AMU)   /* deuteron mass in g    */
#define MASS_T     (3.016042*MASS_AMU)   /* triton   mass in g    */
#define MASS_A     (4.001506*MASS_AMU)   /* alpha    mass in g    */

// DT opacity, from fit to Sesame 
#define OPACITY(T,NE)                                                     \
        ( 0.06*pow((NE)/NA,0.5)*pow(T,-3.5) + 0.14 ) /* cm^-1 */

// Max fractional change to various quantities - used to set time step 
#define MAX_FRAC_CHANGE (0.001)

// Method declarations

struct output_data { 
float *ne; 
float *nd; 
float *nt; 
float *na; 
float *te; 
float *ti; 
float *tr; 
float *time; 
int   *step; 
}; 

void
print_log( const char *fmt, ... );

void 
print_output( struct output_data *od, FILE *fp );

void
e_split_frac(int isplit, float ti, float te, float ne, float *esplit_e); 

/* 
 Burn package. 

 ---------------------------------------
 Physical Units: 

 Time           (sec) 
 Density        (cm^-3)
 Temperature    (keV)
 Reaction rate  (cm^3/s)
 Mass           (g)
 ---------------------------------------
*/ 

int 
main( int argc, char *argv[] ) { 
 float ne, nd, nt, na; 
 float te, ti, tr; 
 float ne_new, nd_new, nt_new, na_new; 
 float te_new, ti_new, tr_new; 
 float *tval, *svbar; 
 float time=0, tstop=1.0e-8, dt; 
 float rmult=1; 
 int   step=0; 
 int   itab, ithick; 
 FILE  *fp; 
 struct output_data odata[1]; 
 int   isplit;

 // Initialize data for diagnostics
 odata->ne   = &ne; 
 odata->nd   = &nd; 
 odata->nt   = &nt; 
 odata->na   = &na; 
 odata->te   = &te; 
 odata->ti   = &ti; 
 odata->tr   = &tr; 
 odata->time = &time; 
 odata->step = &step; 

 // Allocate space for DT reaction rate tables
 ALLOCATE( tval,  NUMDT, float );     
 ALLOCATE( svbar, NUMDT, float );     

 // Read DT reaction rates from file
 fp = fopen( "dt.dat", "r" );
 if ( fp==NULL ) ERROR(( "Couldn't open file dt.dat" )); 
 for ( itab=0; itab<NUMDT; ++itab ) { 
   fscanf( fp, "%e %e", tval+itab, svbar+itab ); 
   // fprintf( stdout, "%e %e\n", tval[itab], svbar[itab] );
 } // for
 fclose( fp ); 


 /////////////////////////////////////////////////////////////////////////// 

 // Physics problem initialization

 // 100 g/cc DT plasma has ne = 2.4 x 10^25 /cc
 // 4.16 g/cc DT plasma has ne = 1.0 x 10^24 /cc

 nd = 0.5e26;
 nt = 0.5e26;
 na = 0; 
 ne = nd + nt + 2.0*na; 
 te = 10.0;
 ti = 10.0;
 tr = 10.0; 
 rmult  = 1.0;     // Reaction rate multiplier
 ithick = 1;       // flag - either 0 or 1
 tstop  = 5.0e-9;  // ne=1.E25/cc
 tstop  = 5.0e-9;  // ne=1.E26/cc
 tstop  = 5.0e-9 ; // ne=1.E24/cc
 isplit = 1;       // Fraley energy splitting: isplit=1
 isplit = 2;       // BPS energy splitting   : isplit=2
 isplit = 3;       // S energy splitting     : isplit=3

 /////////////////////////////////////////////////////////////////////////// 



 // Simulation time advance
 do { 
   double dna_dt, svb;
   float dte_dt_rad=0, dtr_dt_rad=0, kappa, te4_minus_tr4; 
   float dte_dt_fus=0, dti_dt_fus=0, esplit_e; 
   float dte_dt_equ=0, dti_dt_equ=0, nue_ei, nue_de, nue_te, nue_ae; 
   float inv_ni, log_lambda;

#if 1
   // update ion density
   COMPUTE_REACTION_RATE( ti, svb );
   dna_dt = rmult*(svb*nt)*nd;               // parentheses to avoid overflow

   // DEBUG
   //MESSAGE(("dna_dt = %e, svb = %e", dna_dt, svb)); 

   // update ele and rad temperatures from e-rad coupling
   kappa = OPACITY( te, ne );
   te4_minus_tr4 = ( pow(te,4) - pow(tr,4) );
   dte_dt_rad = TWO_THIRDS * CVAC * RAD_CONST * kappa * (-te4_minus_tr4)
                / ( ne * KBOLTZ );
   if ( ithick ) {           // optically thick
     dtr_dt_rad = 0.25 * CVAC * kappa * te4_minus_tr4
                  / ( pow(tr,3) + FLT_MIN );
   } else {                 // optically thin  - Tr = Trad of walls
     dtr_dt_rad = 0; 
   } // if 

   // update ele and ion temperatures from fusion alpha deposition

   e_split_frac( isplit, ti, te, ne, &esplit_e ); // electron splitting factor
   //   fprintf( stdout, "%f \n", esplit_e );     // DEBUG

   dte_dt_fus = TWO_THIRDS * dna_dt * EALPHA * esplit_e 
                / ( ne         * KBOLTZ ); 
   dti_dt_fus = TWO_THIRDS * dna_dt * EALPHA * ( 1.0 - esplit_e ) 
                / ( (nd+nt+na) * KBOLTZ ); 
#endif 

   // electron/ion equilibration
   log_lambda = 5.;                         // FIXME: use real log lambda 
   nue_ei = 1.8e-19 * (   sqrt( MASS_E ) * sqrt( MASS_D ) * 1 * 1 * nd * log_lambda
                          * pow( 1.e3 * (MASS_E*ti + MASS_D*te), -1.5 ) 
                        + sqrt( MASS_E ) * sqrt( MASS_T ) * 1 * 1 * nt * log_lambda
                          * pow( 1.e3 * (MASS_E*ti + MASS_T*te), -1.5 ) 
                        + sqrt( MASS_E ) * sqrt( MASS_A ) * 1 * 4 * na * log_lambda
                          * pow( 1.e3 * (MASS_E*ti + MASS_A*te), -1.5 ) ); 

   nue_de = 1.8e-19 * (   sqrt( MASS_D ) * sqrt( MASS_E ) * 1 * 1 * ne * log_lambda 
                          * pow( 1.e3 * (MASS_D*te + MASS_E*ti), -1.5 ) );

   nue_te = 1.8e-19 * (   sqrt( MASS_T ) * sqrt( MASS_E ) * 1 * 1 * ne * log_lambda 
                          * pow( 1.e3 * (MASS_T*te + MASS_E*ti), -1.5 ) );

   nue_ae = 1.8e-19 * (   sqrt( MASS_A ) * sqrt( MASS_E ) * 4 * 1 * ne * log_lambda 
                          * pow( 1.e3 * (MASS_A*te + MASS_E*ti), -1.5 ) );

   dte_dt_equ = nue_ei * ( ti - te ); 

   // Strange ordering to avoid overflow
   inv_ni = 1.0 / ( nd + nt + na ); 
   dti_dt_equ = ( nue_de*(nd*inv_ni) + nue_te*(nt*inv_ni) + nue_ae*(na*inv_ni) ) * ( te - ti ); 

   // DEBUG
   //printf("nue_de, nue_te, nue_ae = %e %e %e\n", nue_de, nue_te, nue_ae ); 
   //printf("dte_dt_equ, dti_dt_equ = %e %e\n", dte_dt_equ, dti_dt_equ); 


   // Compute increment in various quantities
   BEGIN_PRIMITIVE { 
     int   i; 
     float dnd_dt_tot, dnt_dt_tot, dna_dt_tot, 
           dte_dt_tot, dtr_dt_tot, dti_dt_tot; 
     float dt_max[5] = { FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX }; // Max dt 
                                                                        // change for 
                                                                        // ea. quantity

     dnd_dt_tot = -dna_dt; 
     dnt_dt_tot = -dna_dt; 
     dna_dt_tot = +dna_dt; 
     dte_dt_tot = dte_dt_rad + dte_dt_fus + dte_dt_equ; 
     dtr_dt_tot = dtr_dt_rad; 
     dti_dt_tot =              dti_dt_fus + dti_dt_equ; 

     // Set time step.  
     // Ensure change in nd, nt, te, tr, ti <= MAX_FRAC_CHANGE
     // Use that each quantity is nonnegative-definite. 

#     if 1
     dt_max[0] = MAX_FRAC_CHANGE * nd / ( fabs(dnd_dt_tot) + FLT_MIN );
     dt_max[1] = MAX_FRAC_CHANGE * nt / ( fabs(dnt_dt_tot) + FLT_MIN );
#     endif 

     dt_max[2] = MAX_FRAC_CHANGE * te / ( fabs(dte_dt_tot) + FLT_MIN );
#     if 1
     dt_max[3] = MAX_FRAC_CHANGE * tr / ( fabs(dtr_dt_tot) + FLT_MIN );
#     endif 
     dt_max[4] = MAX_FRAC_CHANGE * ti / ( fabs(dti_dt_tot) + FLT_MIN );
     dt = dt_max[0];
     for ( i=1; i<5; ++i ) {
       if ( dt>dt_max[i] ) {
         dt=dt_max[i];
       } // if
     } // for

     nd   += dnd_dt_tot * dt; 
     nt   += dnt_dt_tot * dt; 
     na   += dna_dt_tot * dt; 
     te   += dte_dt_tot * dt; 
     tr   += dtr_dt_tot * dt; 
     ti   += dti_dt_tot * dt; 
     time += dt; 
     step += 1;
   } END_PRIMITIVE; 

   print_output( odata, stdout ); 

 } while ( time<tstop ); 

 return 0; 
}


void 
print_log( const char *fmt, ... ) {
 va_list ap;
 va_start( ap, fmt );
 vfprintf( stderr, fmt, ap );
 va_end( ap );
 fflush( stderr );
}


void 
print_output( struct output_data *od, FILE *fp ) {
 if ( !fp ) {
   ERROR(( "Bad file pointer passed to print_output()" )); 
 }
 fprintf( fp, "%i %e %e %e %e %e %e %e %e\n", 
          *od->step, *od->time, 
          *od->nd,   *od->nt,   *od->na,   *od->ne, 
          *od->ti,   *od->te,   *od->tr ); 
}



void
e_split_frac(int isplit, float ti, float te, float ne, float *esplit_e_ptr) {
 float esplit_i; 
 float bFR, bS, cne;

 if (isplit == 1) {
   bFR=32.;                      // Fraley offset in keV 
   esplit_i = te / ( te + bFR ); // Fraley energy split
 }
 else if (isplit == 2) {
   cne=-0.005*pow(log10(ne),2) + 0.215*log10(ne) - 0.87;            // BPS fit
   esplit_i=0.5*( 1 + tanh(1.34*( log10(te) - 0.0015*ti - cne )) ); //
 }
 else if (isplit == 3) {
   bS=42.5*3.5/4.;
   esplit_i = te / ( te + bS ); 
 }
 else{
   fprintf( stdout, "%s\n", "isplit=1,2 only" );
   exit(1);
 }

 *esplit_e_ptr = 1.0 - esplit_i;       // electron splitting = 1 - (ion splitting)

 //fprintf( stdout, "%f %f %f\n", *esplit_e_ptr, esplit_i, *esplit_e_ptr + esplit_i ); // DEBUG
}

