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

    class_<CylinderScalingFunctionsWrap<WorldTraits>, bases<CylinderScalingFunctions<WorldTraits> > >("CylinderScalingFunctions", init<PyObject*>() )
        .def("r_right", &CylinderScalingFunctions<WorldTraits>::r_right, &CylinderScalingFunctionsWrap<WorldTraits>::r_right_default)
        ;

    //def( "length_sq", &length_sq<Position> );
    def( "get_dr_dzright_dzleft_to_orthogonal_CylindricalShape", (Position(*)(Position const&, Length const&, Length const&, char*))&get_dr_dzright_dzleft_to_orthogonal_CylindricalShape<WorldTraits> );
}


} // namespace binding
