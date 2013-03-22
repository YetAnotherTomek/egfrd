#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include <boost/python/docstring_options.hpp>
#include "binding_common.hpp"


namespace binding {


void register_cylinder_scaling()
{
    using namespace boost::python;

    docstring_options doc_options;
    doc_options.disable_signatures();

    //def( "length_sq", &length_sq<Position> );
    def( "get_dr_dzright_dzleft_to_orthogonal_CylindricalShape", &get_dr_dzright_dzleft_to_orthogonal_CylindricalShape<WorldTraits> );
}


} // namespace binding
