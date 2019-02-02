/** \file    df_halo.h
    \brief   Distribution functions for the spheroidal component (halo)
    \author  Eugene Vasiliev
    \date    2015-2016
*/
#pragma once
#include "df_base.h"

namespace df{

/// \name   Classes for action-based double power-law distribution function (DF)
///@{

/// Parameters that describe a double power law distribution function.
struct DoublePowerLawParam{
double
    norm,      ///< normalization factor with the dimension of mass
    J0,        ///< break action (defines the transition between inner and outer regions)
    Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
    slopeIn,   ///< power-law index for actions below the break action (Gamma)
    slopeOut,  ///< power-law index for actions above the break action (Beta)
    steepness, ///< steepness of the transition between two asymptotic regimes (eta)
    cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
    coefJrIn,  ///< contribution of radial   action to h(J), controlling anisotropy below J_0 (h_r)
    coefJzIn,  ///< contribution of vertical action to h(J), controlling anisotropy below J_0 (h_z)
    coefJrOut, ///< contribution of radial   action to g(J), controlling anisotropy above J_0 (g_r)
    coefJzOut, ///< contribution of vertical action to g(J), controlling anisotropy above J_0 (g_z)
    rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
    Jphi0;     ///< controls the steepness of rotation and the size of non-rotating core
DoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
    norm(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
    coefJrIn(1), coefJzIn(1), coefJrOut(1), coefJzOut(1), rotFrac(0), Jphi0(0) {}
};

/** General double power-law model.
    The distribution function is given by
    \f$  f(J) = norm / (2\pi J_0)^3
         (1 + (J_0 /h(J))^\eta )^{\Gamma / \eta}
         (1 + (g(J)/ J_0)^\eta )^{-B / \eta }
         \exp[ - (g(J) / J_{cutoff})^\zeta ] \f$,  where
    \f$  g(J) = g_r J_r + g_z J_z + g_\phi |J_\phi|  \f$,
    \f$  h(J) = h_r J_r + h_z J_z + h_\phi |J_\phi|  \f$.
    Gamma is the power-law slope of DF at small J (slopeIn), and Beta -- at large J (slopeOut),
    the transition occurs around J=J0, and its steepness is adjusted by the parameter eta.
    h_r, h_z and h_phi control the anisotropy of the DF at small J (their sum is always taken
    to be unity, so that there are two free parameters -- coefJrIn = h_r, coefJzIn = h_z),
    and g_r, g_z, g_phi do the same for large J (coefJrOut = g_r, coefJzOut = g_z).
    Jcutoff is the threshold of an optional exponential suppression, and zeta measures its strength.
*/
class DoublePowerLaw: public BaseDistributionFunction{
    const DoublePowerLawParam par;  ///< parameters of DF
public:
    /** Create an instance of double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
    DoublePowerLaw(const DoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
    virtual double value(const actions::Actions &J) const;
};

///@}

/// \name   Classes for action-based cored double power-law distribution function (DF)
///@{

/// Parameters that describe a cored double power law distribution function.
struct CoredDoublePowerLawParam : DoublePowerLawParam {
  double
    norm,      ///< normalization factor with the dimension of mass
    J0,        ///< break action (defines the transition between inner and outer regions)
    Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
    slopeIn,   ///< power-law index for actions below the break action (Gamma)
    slopeOut,  ///< power-law index for actions above the break action (Beta)
    steepness, ///< steepness of the transition between two asymptotic regimes (eta)
    cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
    coefJrIn,  ///< contribution of radial   action to h(J), controlling anisotropy below J_0 (h_r)
    coefJzIn,  ///< contribution of vertical action to h(J), controlling anisotropy below J_0 (h_z)
    coefJrOut, ///< contribution of radial   action to g(J), controlling anisotropy above J_0 (g_r)
    coefJzOut, ///< contribution of vertical action to g(J), controlling anisotropy above J_0 (g_z)
    rotFrac,   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
    Jphi0,     ///< controls the steepness of rotation and the size of non-rotating core
    h0,  ///< sets the scale of the almost constant-density core of the final DF
    coreAmp;  ///< temporary for now, should be determined as in Eq. 4 from Cole & Binney 18
CoredDoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
    norm(NAN), J0(NAN), Jcutoff(0), slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
    coefJrIn(1), coefJzIn(1), coefJrOut(1), coefJzOut(1), rotFrac(0), Jphi0(0),
    h0(0), coreAmp(0) {}
};
 
/** Cored double power-law model.
    The distribution function is given by the general double power-law DF
    convolved with a function which moves mass from the inner part to the outside,
    effectively heating the DM core
    coreHeatingScale sets the scale of the almost constant-density core of the final DF
*/
class CoredDoublePowerLaw: public BaseDistributionFunction{
    const CoredDoublePowerLawParam par;  ///< parameters of DF
public:
    /** Create an instance of cored double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
    CoredDoublePowerLaw(const CoredDoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
    virtual double value(const actions::Actions &J) const;
};

///@}
}  // namespace df
