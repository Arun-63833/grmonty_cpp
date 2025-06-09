#pragma once
#include <cmath>  // for M_PI

namespace constants
{
    // Fundamental physical constants (CGS units)
    inline constexpr double EE    = 4.80320680e-10;     // Electron charge
    inline constexpr double CL    = 2.99792458e10;      // Speed of light
    inline constexpr double ME    = 9.1093826e-28;      // Electron mass
    inline constexpr double MP    = 1.67262171e-24;     // Proton mass
    inline constexpr double MN    = 1.67492728e-24;     // Neutron mass
    inline constexpr double AMU   = 1.66053886e-24;     // Atomic mass unit
    inline constexpr double HPL   = 6.6260693e-27;      // Planck constant
    inline constexpr double HBAR  = HPL / (2.0 * M_PI); // Reduced Planck constant
    inline constexpr double KBOL  = 1.3806505e-16;      // Boltzmann constant
    inline constexpr double GNEWT = 6.6742e-8;          // Gravitational constant
    inline constexpr double SIG   = 5.670400e-5;        // Stefan-Boltzmann constant
    inline constexpr double RGAS  = 8.3143e7;           // Ideal gas constant (erg·K⁻¹·mol⁻¹)
    inline constexpr double EV    = 1.60217653e-12;     // Electron volt in erg
    inline constexpr double SIGMA_THOMSON = 0.665245873e-24; // Thomson cross-section
    inline constexpr double JY    = 1.e-23;             // Jansky (flux density unit)

    // Astronomical distance units
    inline constexpr double PC = 3.085678e18;           // Parsec in cm
    inline constexpr double AU = 1.49597870691e13;      // Astronomical unit in cm

    // Time units
    inline constexpr double YEAR = 31536000.0;          // Seconds in a year
    inline constexpr double DAY  = 86400.0;             // Seconds in a day
    inline constexpr double HOUR = 3600.0;              // Seconds in an hour

    // Solar values
    inline constexpr double MSUN = 1.989e33;            // Solar mass
    inline constexpr double RSUN = 6.96e10;             // Solar radius
    inline constexpr double LSUN = 3.827e33;            // Solar luminosity
    inline constexpr double TSUN = 5.78e3;              // Sun's photosphere temperature (K)

    // Earth values
    inline constexpr double MEARTH = 5.976e27;          // Earth mass
    inline constexpr double REARTH = 6.378e8;           // Earth radius

    // Galactic center
    inline constexpr double DSGRA = 8.4e3 * PC;         // Distance to Sgr A*

    // Cosmic background radiation
    inline constexpr double TCBR = 2.726;               // CMB temperature

    // Solar abundances
    inline constexpr double SOLX = 0.70;                // Hydrogen
    inline constexpr double SOLY = 0.28;                // Helium
    inline constexpr double SOLZ = 0.02;                // Metals

    // Math constants
    inline constexpr double FOUR_PI     = 4.0 / M_PI;
    inline constexpr double FOUR_PISQ   = 4.0 / (M_PI * M_PI);
    inline constexpr double THREEPI_TWO = 3.0 * M_PI / 2.0;
}
