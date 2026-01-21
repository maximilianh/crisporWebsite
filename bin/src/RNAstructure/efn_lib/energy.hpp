//
// Created by Max Ward on 27 March 2017.
//

#ifndef EFNMAX_ENERGY_HPP
#define EFNMAX_ENERGY_HPP

namespace efnmax {

/// A value in Kcals/mol.
typedef double kcalmol_t;

/// The increment size of an energy value. Needed since energy is actually an integral value.
/// Currently desgiend to be the same as RNAstructure for compatibility.
const kcalmol_t EnergyIota = 0.1;

/// Standard energy representation. Needs to be an integer type, and stores energy as a multiple of EnergyIota.
typedef int energy_t;

/**
 * @param kcalsmol Energy value is kcals/mol.
 * @return Energy value in increments.
 */
energy_t KCalToEnergy(kcalmol_t kcalsmol);

/**
 * Converts energy increments to kcal/mol.
 * @param e Energy in increments.
 * @return The value in kcal/mol.
 */
kcalmol_t EnergyToKCal(energy_t e);
}

#endif //EFNMAX_ENERGY_HPP
