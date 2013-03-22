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




// TODO the cpp file is not really needed any more
template <typename Ttraits_>
typename Ttraits_::position_type
get_dr_dzright_dzleft_to_CylindricalShape_helper( typename Ttraits_::position_type      const& pos,
                                                  typename Ttraits_::length_type        const& l1,
                                                  typename Ttraits_::length_type        const& l2   )

{
  
    typedef typename Ttraits_::position_type    position_type;
    typedef typename Ttraits_::length_type      length_type;

    position_type return_vector( pos );         // only pass on pos for now
      
    return return_vector;
  
};
                                                  

