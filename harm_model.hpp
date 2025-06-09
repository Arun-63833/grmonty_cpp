#pragma once
#include "decs.hpp" // for struct of_photon, NDIM etc.

namespace harm {

// Instead of `extern`, keep them in a namespace so no macro trick needed
inline double *****econ = nullptr;
inline double *****ecov = nullptr;
inline double ****bcon = nullptr;
inline double ****bcov = nullptr;
inline double ****ucon = nullptr;
inline double ****ucov = nullptr;
inline double ****p    = nullptr;
inline double ***ne    = nullptr;
inline double ***thetae = nullptr;
inline double ***b     = nullptr;

// Function declarations (unchanged for now)
void init_weight_table();
void bl_coord(double *X, double *r, double *th);
void make_zone_centered_tetrads();
void set_units(char *munitstr);
void init_geometry();
void init_harm_data(char *fname);
void init_nint_table();
void init_storage();
double dOmega_func(double x2i, double x2f);

void sample_zone_photon(int i, int j, double dnmax, of_photon *ph);
double interp_scalar(double **var, int i, int j, double coeff[4]);
int get_zone(int *i, int *j, double *dnamx);
void Xtoij(double X[constants::NDIM], int *i, int *j, double del[constants::NDIM]);
void coord(int i, int j, double *X);
void get_fluid_zone(int i, int j, double *Ne, double *Thetae, double *B,
                    double Ucon[constants::NDIM], double Bcon[constants::NDIM]);

} // namespace harm
