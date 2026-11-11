// Solar zenith angle from Gregorian date, elapsed time and geographic position.
// cos_solar_zen returns a double-precision cosine; latitude/longitude are degrees
// north/east and current_time is seconds since the supplied start date.
// References: Clayson and Curry (1996),JGR101(C12),28515–28528;
// Liou, An Introduction to Atmospheric Radiation,table2.2,equation2.2.10.
#pragma once
#include "clubb_real.hpp"

namespace clubb_kokkos {

double cos_solar_zen(int day, int month, int year, double current_time, // In
                     real_t lat_in_degrees, real_t lon_in_degrees       // In
);

} // namespace clubb_kokkos
