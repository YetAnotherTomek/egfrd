#ifndef BINDING_CYLINDRICAL_SURFACE_HPP
#define BINDING_CYLINDRICAL_SURFACE_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"
#include "Defs.hpp"

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_cylindrical_surface_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(
            name, init<typename impl_type::structure_name_type,
                       typename impl_type::structure_type_id_type,
                       typename impl_type::structure_id_type,
                       typename impl_type::shape_type>())
        .def(init<typename impl_type::structure_name_type,
                       typename impl_type::structure_type_id_type,
                       typename impl_type::structure_id_type,
                       typename impl_type::shape_type, Real, Real>())
        .add_property("shape",
            make_function((typename impl_type::shape_type const&(impl_type::*)()const)&impl_type::shape, return_value_policy<reference_existing_object>()))
        .add_property("growth_rate",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    Real,
                    &impl_type::growth_rate,
                    &impl_type::growth_rate>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    impl_type,
                    Real,
                    &impl_type::growth_rate,
                    &impl_type::growth_rate>::set))
        ;
}

} // namespace binding

#endif /* BINDING_CYLINDRICAL_SURFACE_HPP */
