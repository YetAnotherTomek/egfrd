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




/*** LOGGER ***/        // FIXME
//Logger& loclog_(Logger::get_logger("ecell.CylinderScaling"));


// Declare a base class with bogus scaling methods which will
// be overridden in Python via a derived class
template <typename Ttraits_>
class CylinderScalingFunctions
{
 public:
   
    typedef typename Ttraits_::length_type      length_type;

    // Here we define virtual default functions which can be overridden by both C++ and Python
    // classes derived from this base class
    virtual length_type r_right(length_type z) { return 0.0; };
    virtual length_type z_right(length_type r) { return 0.0; };
    virtual length_type r_left(length_type z)  { return 0.0; };
    virtual length_type z_left(length_type r)  { return 0.0; };
    virtual ~CylinderScalingFunctions() {};
};


// These are testing methods to call methods of the wrapped class
template <typename Ttraits_>
typename Ttraits_::length_type calls_r_right(CylinderScalingFunctions<Ttraits_> *CSF, typename Ttraits_::length_type z) { return CSF->r_right(z); };

template <typename Ttraits_>
typename Ttraits_::length_type calls_z_right(CylinderScalingFunctions<Ttraits_> *CSF, typename Ttraits_::length_type r) { return CSF->z_right(r); };

template <typename Ttraits_>
typename Ttraits_::length_type calls_r_left(CylinderScalingFunctions<Ttraits_> *CSF, typename Ttraits_::length_type z) { return CSF->r_left(z); };

template <typename Ttraits_>
typename Ttraits_::length_type calls_z_left(CylinderScalingFunctions<Ttraits_> *CSF, typename Ttraits_::length_type r) { return CSF->z_left(r); };

// Wrapping Code

// Define a dispatcher class derived from the base class which takes
// a Python class as a parameter and overrides the methods of the
// base class with the Python methods via boost::python::call_method
template <typename Ttraits_>
class CylinderScalingFunctionsWrap : public CylinderScalingFunctions<Ttraits_>
{
  
  typedef typename Ttraits_::length_type      length_type;
  
  public:       

    // The constructor takes a (potentially) derived Python version of
    // class CylinderScalingFunctions as an argument
    CylinderScalingFunctionsWrap(PyObject* self_) : self(self_) {};    
    ~CylinderScalingFunctionsWrap() {};

    // Default implementation, for when functions are not overridden
    length_type r_right_default(length_type z) { return this->CylinderScalingFunctions<Ttraits_>::r_right(z); };
    length_type z_right_default(length_type r) { return this->CylinderScalingFunctions<Ttraits_>::z_right(r); };
    length_type r_left_default(length_type z)  { return this->CylinderScalingFunctions<Ttraits_>::r_left(z); };
    length_type z_left_default(length_type r)  { return this->CylinderScalingFunctions<Ttraits_>::z_left(r); };
    // Dispatch implementation, using boost::python::call_method on the PyObject* passed to the constructor
    length_type r_right(length_type z) { return boost::python::call_method<length_type>(self, "r_right", z); };
    length_type z_right(length_type r) { return boost::python::call_method<length_type>(self, "z_right", r); };
    length_type r_left(length_type z)  { return boost::python::call_method<length_type>(self, "r_left", z); };
    length_type z_left(length_type r)  { return boost::python::call_method<length_type>(self, "z_left", r); };
    
  public:
    
    // Store a pointer to the Python object in property this->self
    PyObject* self;
};





template <typename Ttraits_>
class CylinderScalingHelperTools
{
 
  typedef typename Ttraits_::length_type      length_type;
  typedef typename Ttraits_::length_type      position_type;
  
  typedef int                                 direction_type;
  
  // define scaling helper function type to later define an array
  // of scaling functions with direction_type index
  typedef length_type (CylinderScalingHelperTools<Ttraits_>::* scaling_function_type)(length_type);
  
  
  public:       // methods
    
    CylinderScalingHelperTools(CylinderScalingFunctionsWrap<Ttraits_> *CSF_, Real scale_angle_) : 
      CSF(CSF_), scale_angle(scale_angle_)
    {
         // Assign the scaling functions to the array for the two different directions
         // direction = -1 (direction_index = 0)
         r1_function[0] = &CylinderScalingHelperTools<Ttraits_>::r_left;
         z1_function[0] = &CylinderScalingHelperTools<Ttraits_>::z_left;
         z2_function[0] = &CylinderScalingHelperTools<Ttraits_>::z_right;
         // direction = +1 (direction_index = 1)
         r1_function[1] = &CylinderScalingHelperTools<Ttraits_>::r_right;
         z1_function[1] = &CylinderScalingHelperTools<Ttraits_>::z_right;
         z2_function[1] = &CylinderScalingHelperTools<Ttraits_>::z_left;
         
         tan_scale_angle = std::tan(scale_angle);
    };
    
    ~CylinderScalingHelperTools() {};

    // master method that calls all the others in the right order
    //inline static position_type get_dr_dzright_dzleft_to_CylindricalShape();
    // helper methods
    //inline static void determine_direction_and_coordinate_system();
    inline static position_type get_dr_dzright_dzleft_to_orthogonal_CylindricalShape();
    inline static position_type get_dr_dzright_dzleft_to_parallel_CylindricalShape();    

    inline static position_type get_dr_dzright_dzleft_to_CylindricalShape()
    {
        return position_type();
    };
    
    length_type test_r1_function(length_type z)  // TESTING
    {
        this->direction_index = int(scale_angle == 1.0 ? 1 : 0);
        
        return (this->*r1_function[this->direction_index])(z);
    };
    
    length_type test_z1_function(length_type r)  // TESTING
    {
        this->direction_index = int(scale_angle == 1.0 ? 1 : 0);
        
        return (this->*z1_function[this->direction_index])(r);
    };
    
  private:

    inline static void determine_direction_and_coordinate_system()
    {
    };
    
    // Wrappers for cylinder scaling functions passed as member functions of CSF
    // This overhead is necessary because C++ cannot directly form a pointer to  
    // the bound member functions of CSF
    length_type r_left (length_type z){  return CSF->r_left(z);  };
    length_type r_right(length_type z){  return CSF->r_right(z); };
    length_type z_left (length_type r){  return CSF->z_left(r);  };
    length_type z_right(length_type r){  return CSF->z_right(r); };
    
  public:       // properties
    
    CylinderScalingFunctionsWrap<Ttraits_> *CSF;
    
    direction_type      direction;
    int                 direction_index; // to address the methods in the scaling function pointer arrays
                                         // has to start from zero, so we have to map direction=-1 to direction_index=0
    
    position_type       orientation_vector;
    position_type       ref_to_shell_vector;    
    length_type         ref_to_shell_z;
    
    Real                scale_angle;
    Real                tan_scale_angle;
    
    // array of scaling methods - we have to call different ones depending on the direction
    // determined initially
    scaling_function_type r1_function[2];
    scaling_function_type z1_function[2];
    scaling_function_type z2_function[2];

};



// OLDER STUFF; TODO Remove when the rest works
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
    //double result;
    
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
