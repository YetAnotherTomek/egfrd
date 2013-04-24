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

    // Bind the overridden CylinderScalingFunctions class into Python; the derived class takes a Python class as an argument
    // and overrides the bogus methods defined in the base class by the Python methods
    class_<CylinderScalingFunctions<WorldTraits> >("CylinderScalingFunctionsBase");
    
    class_<CylinderScalingFunctionsWrap<WorldTraits>, bases<CylinderScalingFunctions<WorldTraits> >, boost::shared_ptr<CylinderScalingFunctionsWrap<WorldTraits> >, boost::noncopyable >
        ("CylinderScalingFunctions", init<PyObject*>() )   // here we also specify that the constructor gets a pointer to the (possibly derived) Python version of the wrapped class
        .def("r_right", &CylinderScalingFunctions<WorldTraits>::r_right, &CylinderScalingFunctionsWrap<WorldTraits>::r_right_default)
        .def("z_right", &CylinderScalingFunctions<WorldTraits>::z_right, &CylinderScalingFunctionsWrap<WorldTraits>::z_right_default)
        .def("r_left", &CylinderScalingFunctions<WorldTraits>::r_left, &CylinderScalingFunctionsWrap<WorldTraits>::r_left_default)
        .def("z_left", &CylinderScalingFunctions<WorldTraits>::z_left, &CylinderScalingFunctionsWrap<WorldTraits>::z_left_default)
        ;
        
    // Bind the functions that call the class methods after the double-wrapping magic
    def("calls_r_right", &calls_r_right<WorldTraits>);
    def("calls_z_right", &calls_z_right<WorldTraits>);
    def("calls_r_left", &calls_r_left<WorldTraits>);
    def("calls_z_left", &calls_z_left<WorldTraits>);
    
    class_<CylinderScalingHelperTools<WorldTraits>, boost::shared_ptr<CylinderScalingHelperTools<WorldTraits> >, boost::noncopyable >
        ("CylinderScalingHelperTools", 
          init<CylinderScalingFunctionsWrap<WorldTraits>*, Position, Position, Position, Position, Position, Position, int, Position, Position>() 
        )
        .def("get_dr_dzright_dzleft_to_CylindricalShape", &CylinderScalingHelperTools<WorldTraits>::get_dr_dzright_dzleft_to_CylindricalShape)
        .def("test_z1_function", &CylinderScalingHelperTools<WorldTraits>::test_z1_function)
        .def("test_r1_function", &CylinderScalingHelperTools<WorldTraits>::test_r1_function)
        ;    
}


} // namespace binding

