#ifndef CYLINDERSCALING_HPP
#define CYLINDERSCALING_HPP

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include <math.h>
//#include <Python.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/call_method.hpp>
#include <boost/utility.hpp>


#include "findRoot.hpp"
#include "funcSum.hpp"
#include "geometry.hpp"

#include "Logger.hpp"



/*** LOGGER ***/
//Logger& loclog_(Logger::get_logger("ecell.CylinderScaling"));



template <typename Ttraits_>
class CylinderScalingFunctions
{
 public:
   
    typedef typename Ttraits_::length_type      length_type;
    
    virtual length_type r_right(length_type z) { return 0.0; };
    virtual length_type z_right(length_type r) { return 0.0; };
    virtual length_type r_left(length_type z)  { return 0.0; };
    virtual length_type z_left(length_type r)  { return 0.0; };
    //virtual ~CylinderScalingFunctions();      // FIXME why does it not link properly with this included? -> Check for memory leaks!
};


template <typename Ttraits_>
typename Ttraits_::length_type calls_r_right(CylinderScalingFunctions<Ttraits_> const& CSF, typename Ttraits_::length_type z) { return CSF.r_right(z); };

//
// Wrapping Code
//

// Dispatcher class
template <typename Ttraits_>
class CylinderScalingFunctionsWrap : public CylinderScalingFunctions<Ttraits_>
{
  public:
    
    typedef typename Ttraits_::length_type      length_type;
    
    // Store a pointer to the Python object
    CylinderScalingFunctionsWrap(PyObject* self_) : self(self_) {};
    PyObject* self;

    // Default implementation, for when functions are not overridden
    length_type r_right_default(length_type z) { return this->CylinderScalingFunctions<Ttraits_>::r_right(z); };
    length_type z_right_default(length_type r) { return this->CylinderScalingFunctions<Ttraits_>::z_right(r); };
    length_type r_left_default(length_type z)  { return this->CylinderScalingFunctions<Ttraits_>::r_left(z); };
    length_type z_left_default(length_type r)  { return this->CylinderScalingFunctions<Ttraits_>::z_left(r); };
    // Dispatch implementation
    length_type r_right(length_type z) { return boost::python::call_method<length_type>(self, "r_right", z); };
    length_type z_right(length_type r) { return boost::python::call_method<length_type>(self, "z_right", r); };
    length_type r_left(length_type z)  { return boost::python::call_method<length_type>(self, "r_left", z); };
    length_type z_left(length_type r)  { return boost::python::call_method<length_type>(self, "z_left", r); };
};







template <typename Ttraits_>
static typename Ttraits_::position_type
get_dr_dzright_dzleft_to_orthogonal_CylindricalShape( typename Ttraits_::position_type      const& pos,
                                                      typename Ttraits_::length_type        const& l1,
                                                      typename Ttraits_::length_type        const& l2,   
                                                      char*                                        f1 )
{
  
    typedef typename Ttraits_::position_type    position_type;
    typedef typename Ttraits_::length_type      length_type;
    
    length_type l3(0.0);
    double result;
    
//     std::cout << "Arguments: pos=" << pos 
//               << " l1=" << l1
//               << " l2=" << l2
//               << " f1=" << f1
//               << std::endl;
//     
//     PyObject *pModuleName, *pModule, *pFunc, *pValue, *pArgs, *pResult;
// 
//     // Initialize the Python Interpreter
//     Py_Initialize();
// 
//     pModuleName = PyString_FromString("shells");    
//     //pFuncName   = PyString_FromString(f1);
//     pModule = PyImport_Import(pModuleName);
//     if(pModule) std::cout << "Imported Python module." << std::endl;
//     // pDict is a borrowed reference 
//     //pDict = PyModule_GetDict(pModule);
//     // pFunc is also a borrowed reference 
//     std::cout << "Trying to import function " << f1 << std::endl;
//     //pFunc = PyDict_GetItemString(pDict, f1);
//     pFunc = PyObject_GetAttrString(pModule, f1);
//     std::cout << "Imported Python function." << std::endl;       
//     
//     if (pFunc && PyCallable_Check(pFunc)) 
//     {
//         std::cout << "Building arguments for " << f1 << std::endl;
//         //pArgs = PyTuple_New(1);
//         pArgs = PyTuple_Pack(1, PyFloat_FromDouble(l1));
//         pValue = PyFloat_FromDouble((double)l1);
//         double reconverted( PyFloat_AsDouble(pValue) );
//         std::cout << "Converted argument evaluates to " << reconverted << std::cout;
//         //PyTuple_SetItem(pArgs, 0, pValue);
// 
//         std::cout << "Trying to call function " << f1 << std::endl;
//         pResult = PyObject_CallObject(pFunc, pArgs);
//         
//         if (pResult != NULL) {
//             std::cout << "Converting " << pResult << std::endl;
//             result = PyFloat_AsDouble(pResult);
//             l3 = (length_type)PyFloat_AsDouble(pResult);
//             std::cout << "Converted result = " << result << ", l3 = " << l3 << std::endl;
//         }
//     }
//     else 
//     {
//         std::cout << "Cannot execute Python function." << std::endl;
//         PyErr_Print();
//     }


    std::cout << "Constructing result" << std::endl;
    const position_type return_vector( multiply(pos, l3) );         // TESTING some test output to check functionality of geometry and type conversion
    std::cout << "Result is " << return_vector << std::endl;
    
//     // Clean up    
//     Py_DECREF(pModule);
//     Py_DECREF(pModuleName);
//     Py_DECREF(pFunc);
//     Py_DECREF(pValue);
//     Py_DECREF(pArgs);    
//     Py_DECREF(pResult);
//     
//         
//     // Finish the Python Interpreter
//     std::cout << "Shutting down python interpereter" << std::endl;    
//     Py_Finalize();


    std::cout << "Returning" << std::endl;
    return return_vector;
  
};














#endif /* CYLINDERSCALING_HPP */
