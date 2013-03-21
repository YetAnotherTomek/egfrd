#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include <math.h>

#include "findRoot.hpp"
#include "CylinderScaling.hpp"





template <typename Ttraits_>
std::pair<typename Ttraits_::length_type, std::pair<typename Ttraits_::length_type, typename Ttraits_::length_type> >
get_dr_dzright_dzleft_to_CylindricalShape_helper( typename Ttraits_::position_type      const& pos,
                                                  typename Ttraits_::length_type        const& l1,
                                                  typename Ttraits_::length_type        const& l2   )
{

    // TESTING
    return std::make_pair( l1*l2, return std::make_pair(l1, l2) );
  
};
                                                  

