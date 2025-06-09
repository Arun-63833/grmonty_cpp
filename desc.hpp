
#include<iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include "constants.hpp"


namespace constants {
    inline constexpr int NDIM = 4; // dimensions in space (3 space + time)
    inline constexpr int NPRIM = 8; // primitive variables like density, velocity etc

    //maximum and minimum frequency of super photons
    inline constexpr double NUMIN = 1.e9;
    inline constexpr double NUMAX = 1.e16;

    //maximum and minimum electron temp in the units of rest mass energy
    inline constexpr double THETAE_MAX = 1000.;
    inline constexpr double THETAE_MIN = 0.3;
    
    
    inline constexpr double TP_OVER_TE = 3.0; // ion-to-electron temperature ratio
    inline constexpr double WEIGHT_MIN = 1.e31; // smallest weight of proton considered for tracking 

    inline constexpr int KRHO = 0;
    //velocity 1 , 2 , 3
    inline constexpr int UU   = 1;
    inline constexpr int U1   = 2;
    inline constexpr int U2   = 3;
    inline constexpr int U3   = 4;

    //magnetic field 1 , 2 , 3
    inline constexpr int B1   = 5;
    inline constexpr int B2   = 6;
    inline constexpr int B3   = 7;

    inline constexpr double SMALL = 1.e-40;
    inline constexpr double MMW = 0.5;

    inline constexpr int N_ESAMP = 200;
    inline constexpr int N_EBINS = 200;
    inline constexpr int N_THBINS = 6;
}

struct of_photon {
	double X[constants::NDIM]; // position in space time
	double K[constants::NDIM]; // energy and direction of photon
	double dKdlam[constants::NDIM]; //Derivative of momentum w.r.t affine parameter
	double w; //number of physical proton rep by the super proton
	
    //quantities tracked during photon propagation
    double E;
	double L;
	double X1i;
	double X2i;
	double tau_abs;
	double tau_scatt;
	double ne0;
	double thetae0;
	double b0;
	double E0;
	double E0s;

    //number of scatterings
	int nscatt;
};

struct of_geom {
	double gcon[constants::NDIM][constants::NDIM]; // contravariant metric tensor gⁱʲ
	double gcov[constants::NDIM][constants::NDIM]; // covariant metric tensor gᵢⱼ
	double g; // determinant of metric
};

struct of_spectrum {
	double dNdlE;
	double dEdlE;
	double nph;
	double nscatt;
	double X1iav;
	double X2isq;
	double X3fsq;
	double tau_abs;
	double tau_scatt;
	double ne0;
	double thetae0;
	double b0;
	double E0;
};

struct of_grid {
	struct of_spectrum spec[constants::NDIM];
	double th, phi;
	int nlist;
	int *in;
};


/** global variables **/
/** model independent */
extern gsl_rng *r;

extern double F[constants::NDIM + 1], wgt[constants::NDIM + 1];

extern int Ns;
extern int N_superph_recorded, N_scatt;

/* HARM model globals */
extern struct of_geom **geom;
extern int N1, N2, N3;
extern int n_within_horizon;


/* some coordinate parameters */
extern double a;
extern double R0, Rin, Rh, Rout, Rms;
extern double hslope;
extern double startx[constants::NDIM], stopx[constants::NDIM], dx[constants::NDIM];
extern double dlE, lE0;
extern double gam;
extern double dMsim;

extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;

extern double max_tau_scatt, Ladv, dMact, bias_norm;

/* some useful macros */
#define DLOOP  for(k=0;k< constants::NDIM ;k++)for(l=0;l< constants::NDIM ;l++)
#define INDEX(i,j,k)	(constants::NPRIM*( (k) + N3*((j) + N2*(i))))
#define MYSIN(x,sx) 	{							\
			double _xp = (x)-M_PI; 					\
			double _yp = _xp*(FOUR_PI - FOUR_PISQ*fabs(_xp)); 	\
			sx = -_yp*(0.225*fabs(_yp)+0.775);			\
			}
#define MYCOS(x,cx) 	{							\
			double _xp = (x)-THREEPI_TWO; 					\
			_xp += (_xp<-M_PI)*2.*M_PI; 				\
			double _yp = _xp*(FOUR_PI - FOUR_PISQ*fabs(_xp));		\
			cx = _yp*(0.225*fabs(_yp)+0.775);			\
			}

/** model-independent subroutines **/
/* core monte carlo/radiative transport routines */
void track_super_photon(struct of_photon *ph);
void record_super_photon(struct of_photon *ph);
void report_spectrum(int N_superph_made);
void scatter_super_photon(struct of_photon *ph, struct of_photon *php,
			  double Ne, double Thetae, double B,
			  double Ucon[constants::NDIM], double Bcon[constants::NDIM],
			  double Gcov[constants::NDIM][constants::NDIM]);

/* OpenMP specific functions */
void omp_reduce_spect(void);

/* MC/RT utilities */
void init_monty_rand(int seed);
double monty_rand(void);

/* geodesic integration */
void init_dKdlam(double X[], double Kcon[], double dK[]);
void push_photon_ham(double X[constants::NDIM], double Kcon[][constants::NDIM], double dl[]);
void push_photon(double X[constants::NDIM], double Kcon[constants::NDIM], double dKcon[constants::NDIM],
		 double dl, double *E0, int n);
void push_photon4(double X[constants::NDIM], double Kcon[constants::NDIM], double dKcon[constants::NDIM],
		  double dl);
void push_photon_cart(double X[constants::NDIM], double Kcon[constants::NDIM],
		      double dKcon[constants::NDIM], double dl);
double stepsize(double X[constants::NDIM], double K[constants::NDIM]);
void push_photon_gsl(double X[constants::NDIM], double Kcon[constants::NDIM], double dl);
int geodesic_deriv(double t, const double y[], double dy[], void *params);
void interpolate_geodesic(double Xi[], double X[], double Ki[], double K[],
			  double frac, double del_l);

/* basic coordinate functions supplied by grmonty */
void boost(double k[constants::NDIM], double p[constants::NDIM], double ke[constants::NDIM]);
void lower(double *ucon, double Gcov[constants::NDIM][constants::NDIM], double *ucov);
double gdet_func(double gcov[][constants::NDIM]);  /* calculated numerically */
void coordinate_to_tetrad(double Ecov[constants::NDIM][constants::NDIM], double K[constants::NDIM],
			  double K_tetrad[constants::NDIM]);
void tetrad_to_coordinate(double Ecov[constants::NDIM][constants::NDIM], double K_tetrad[constants::NDIM],
			  double K[constants::NDIM]);
double delta(int i, int j);
void normalize(double Ucon[constants::NDIM], double Gcov[constants::NDIM][constants::NDIM]);
void normalize_null(double Gcov[constants::NDIM][constants::NDIM], double K[constants::NDIM]);
void make_tetrad(double Ucon[constants::NDIM], double Bhatcon[constants::NDIM],
		 double Gcov[constants::NDIM][constants::NDIM], double Econ[constants::NDIM][constants::NDIM],
		 double Ecov[constants::NDIM][constants::NDIM]);

/* functions related to basic radiation functions & physics */
	/* physics-independent */
double get_fluid_nu(double X[4], double K[4], double Ucov[constants::NDIM]);
double get_bk_angle(double X[constants::NDIM], double K[constants::NDIM], double Ucov[constants::NDIM],
		    double Bcov[constants::NDIM], double B);
double alpha_inv_scatt(double nu, double thetae, double Ne);
double alpha_inv_abs(double nu, double thetae, double Ne, double B,
		     double theta);
double Bnu_inv(double nu, double thetae);
double jnu_inv(double nu, double thetae, double ne, double B,
	       double theta);

	/* thermal synchrotron */
double jnu_synch(double nu, double Ne, double Thetae, double B,
		 double theta);
double int_jnu(double Ne, double Thetae, double Bmag, double nu);
void init_emiss_tables(void);
double F_eval(double Thetae, double Bmag, double nu);
double K2_eval(double Thetae);

	/* compton scattering */
void init_hotcross(void);
double total_compton_cross_lkup(double nu, double theta);
double klein_nishina(double a, double ap);
double kappa_es(double nu, double theta);
void sample_electron_distr_p(double k[constants::NDIM], double p[constants::NDIM], double theta);
void sample_beta_distr(double theta, double *gamma_e, double *beta_e);
double sample_klein_nishina(double k0);
double sample_thomson(void);
double sample_mu_distr(double beta_e);
double sample_y_distr(double theta);
void sample_scattered_photon(double k[constants::NDIM], double p[constants::NDIM],
			     double kp[constants::NDIM]);

/** model dependent functions required by code: these 
   basic interfaces define the model **/

/* physics related */
void init_model(char *args[]);
void make_super_photon(struct of_photon *ph, int *quit_flag);
double bias_func(double Te, double w);
void get_fluid_params(double X[constants::NDIM], double gcov[constants::NDIM][constants::NDIM], double *Ne,
		      double *Thetae, double *B, double Ucon[constants::NDIM],
		      double Ucov[constants::NDIM], double Bcon[constants::NDIM],
		      double Bcov[constants::NDIM]);
int stop_criterion(struct of_photon *ph);
int record_criterion(struct of_photon *ph);

/* coordinate related */
void get_connection(double *X, double lconn[][constants::NDIM][constants::NDIM]);
void gcov_func(double *X, double gcov[][constants::NDIM]);
void gcon_func(double *X, double gcon[][constants::NDIM]);

