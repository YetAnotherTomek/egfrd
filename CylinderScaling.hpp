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
#include <Python.h>

#include "findRoot.hpp"
#include "funcSum.hpp"
#include "geometry.hpp"

#include "Logger.hpp"



/*** LOGGER ***/
//Logger& loclog_(Logger::get_logger("ecell.CylinderScaling"));



template <typename Ttraits_>
typename Ttraits_::position_type
get_dr_dzright_dzleft_to_orthogonal_CylindricalShape( typename Ttraits_::position_type      const& pos,
                                                      typename Ttraits_::length_type        const& l1,
                                                      typename Ttraits_::length_type        const& l2,   
                                                      char*                                        f1 )
{
  
    typedef typename Ttraits_::position_type    position_type;
    typedef typename Ttraits_::length_type      length_type;
    
    length_type l3(0.0);
    
    PyObject *pModuleName, *pModule, *pDict, *pFuncName, *pFunc, *pValue;

    // Initialize the Python Interpreter
    Py_Initialize();

    pModuleName = PyString_FromString("shells");    
    //pFuncName   = PyString_FromString(f1);
    pModule = PyImport_Import(pModuleName);
    if (pModule) std::cout << "Imported Python module." << std::endl;
    // pDict is a borrowed reference 
    //pDict = PyModule_GetDict(pModule);
    // pFunc is also a borrowed reference 
    std::cout << "Trying to import function " << f1 << std::endl;
    //pFunc = PyDict_GetItemString(pDict, f1);
    pFunc = PyObject_GetAttrString(pModule, f1);
    std::cout << "Imported Python function." << std::endl;       
    
    if (pFunc && PyCallable_Check(pFunc)) 
    {
        std::cout << "Building arguments for " << f1 << std::endl;
        PyObject* pArgs = PyTuple_Pack(1,PyFloat_FromDouble(l1));
        std::cout << "Trying to call function " << f1 << std::endl;
        PyObject* pResult = PyObject_CallObject(pFunc, pArgs);
        
        std::cout << "Converting " << pResult << std::endl;
        l3 = PyFloat_AsDouble(pResult);
    }
    else 
    {
        std::cout << "Cannot execute Python function." << std::endl;
        PyErr_Print();
    }


    std::cout << "Constructing result" << std::endl;
    position_type return_vector( multiply( pos, l3) );         // TESTING some test output to check functionality of geometry and type conversion

    
    // Clean up
//     Py_DECREF(pModule);
//     Py_DECREF(pModuleName);
    // Finish the Python Interpreter
    Py_Finalize();


    
    return return_vector;
  
};


















#endif /* CYLINDERSCALING_HPP */