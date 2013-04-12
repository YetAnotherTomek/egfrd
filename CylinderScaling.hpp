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
#include "exceptions.hpp"

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
  typedef typename Ttraits_::position_type    position_type;  
  typedef int                                 direction_type;
  typedef Real                                angle_type; // TODO could be improved
  
  enum collision_type { BARREL_HITS_FLAT, EDGE_HITS_EDGE, BARREL_HITS_EDGE, FLAT_HITS_BARREL, EDGE_HITS_BARREL, BARREL_HITS_BARREL };
  
  // define scaling helper function type to later define an array
  // of scaling functions with direction_type index
  typedef length_type (CylinderScalingHelperTools<Ttraits_>::* scaling_function_pt_type)(length_type);
  
  
  public:       // METHODS          
        
    inline position_type get_dr_dzright_dzleft_to_CylindricalShape()
    {
        // The master method that calls downstream helper methods in the right order
        // It returns the new testShell dimensions after scaling against otherShell,
        // after first determining whether the two shells are oriented in parallel 
        // or orthogonally
        
        length_type relative_orientation( dot_product(testShell_orientation_vector, otherShell_orientation_vector) );
        
        if( feq(relative_orientation, 1.0, 1e-3 ) )
          return get_dr_dzright_dzleft_to_parallel_CylindricalShape();
        
        else if( feq(relative_orientation, 0.0, 1e-3 ) )
          return get_dr_dzright_dzleft_to_orthogonal_CylindricalShape();
        
        else
          throw unsupported("Shells are neither parallel nor orthogonal in CylindricaltestShell scaling routine.");
        
    };
    
    // TESTING methods; to test whether scaling functions passed to this class
    // via CylinderScalingFunctionsWrap class are correctly invoked from Python
    length_type test_r1_function(length_type z)  // TESTING
    {        
        return (this->*r1_function[this->di])(z);
    };
    
    length_type test_z1_function(length_type r)  // TESTING
    {        
        return (this->*z1_function[this->di])(r);
    };
    
    
  private:          // HELPER METHODS

    inline void determine_direction()
    {
        // For now this is still done in Python
    };
    
    // Collision of orthogonal cylinders
    inline position_type get_dr_dzright_dzleft_to_orthogonal_CylindricalShape()
    {
        construct_local_coordinate_system();
        
        return get_dr_dzright_dzleft_for_specific_collision( determine_collision_type() );
    };
    
    // Helper method to construct the local coordinate system for orthogonal scaling
    inline void construct_local_coordinate_system()
    {
    };
    
    // Helper method to determine how precisely the testShell hits the otherShell
    inline collision_type determine_collision_type()
    {
    };
    
    // Helper method that calculates the new lengths for the specific collision determined above
    inline position_type get_dr_dzright_dzleft_for_specific_collision(collision_type collision)
    {
    };
    
    // Collision of parallel cylinders
    inline position_type get_dr_dzright_dzleft_to_parallel_CylindricalShape()
    {
        // Calculates the new testShell dimensions when scaled with respect
        // to another shell that is parallel in orientation
        
        length_type r_new(r), z1_new(z1), z2_new(z2);
        
        // calculate ref_to_shell_r/z in the cylindrical coordinate system on the right/left side
        position_type ref_to_shell_z_vec( multiply(testShell_orientation_vector, ref_to_shell_z) );
        position_type ref_to_shell_r_vec( subtract(ref_to_shell_vec, ref_to_shell_z_vec) );
        length_type   ref_to_shell_r( length(ref_to_shell_r_vec) );     // the radius is always positive
        length_type   ref_to_shell_z( ref_to_shell_z * direction );     // ref_to_shell_z is positive 
                                                                        // on the scaling side (right/left)

        // calculate the distances in r/z from the scaling center to the shell
        length_type scale_center_to_shell_r( ref_to_shell_r - scale_center_r );
        length_type scale_center_to_shell_z( ref_to_shell_z - scale_center_z );

        // get angles
        angle_type to_edge_angle( std::atan( (scale_center_to_shell_r - otherShell_radius) / 
                                                 (scale_center_to_shell_z - otherShell_hl)
                                           ) );

        if(scale_center_to_shell_z - otherShell_hl < 0.0 )
            to_edge_angle += M_PI;      // if the shell was too much to the side we correct the angle to be positive
        // elif: a negative angle was caused by a negative scale_center_to_shell we want a negative angle -> do nothing
        // otherwise: shell_angle is ok -> do nothing

        if(to_edge_angle <= scale_angle)
        {   // shell hits the scaling cylinder on top
            z1_new = std::min(z1, (ref_to_shell_z - otherShell_hl) );
            r_new  = std::min(r,  (this->*r1_function[this->di])(z1_new) );
                     // TODO if z1 hasn't changed we also don't have to recalculate this
        }
        else
        {   // shell hits the scaling cylinder on the radial side
            r_new  = std::min(r, (ref_to_shell_r - otherShell_radius) );
            z1_new = std::min(z1, (this->*z1_function[this->di])(r_new) );
        }
        
        return create_vector<position_type>(r_new, z1_new, z2_new);
    };
    
    
    // Wrappers for cylinder scaling functions passed as member functions of CSF
    // This overhead is necessary because C++ cannot directly form a pointer to  
    // the bound member functions of CSF
    length_type r_left (length_type z){  return CSF->r_left(z);  };
    length_type r_right(length_type z){  return CSF->r_right(z); };
    length_type z_left (length_type r){  return CSF->z_left(r);  };
    length_type z_right(length_type r){  return CSF->z_right(r); };
    
    
  public:       // PROPERTIES
    
    CylinderScalingFunctionsWrap<Ttraits_> *CSF;        
    
    // We pack some of the info into position vectors. This is to limit the no. of arguments.
    position_type       testShell_reference_point;
    position_type       testShell_orientation_vector;
    position_type       testShell_dimensions;   // first entry radius r, second z1, third z2
    position_type       otherShell_position_t;
    position_type       otherShell_orientation_vector;
    position_type       otherShell_dimensions;  // 1st entry radius r, 2nd half_length, 3r http://www.imaginefilmfestival.nl/d unused
    position_type       scale_center_info;      // 1st entry r-coord. of the scale center, 2nd z-coord.
    position_type       scale_angle_info;       // 1st entry scale_angle, 2nd entry tan(scale_angle), 3rd unused;
                                                // tan(scale_angle) is calculated in Python already at shell 
                                                // construction because tan() is quite expensive, so passing on 
                                                // is better than recalculating.
    // testShell dimensions
    length_type         r;  // radius
    length_type         z1; // the height that is scaled
    length_type         z2; // the height on the opposite side of the shell (unscaled)
    // otherShell dimensions
    length_type         otherShell_radius;
    length_type         otherShell_hl;  // half-length
    
    // Scaling info follows
    direction_type      direction;
    int                 di; // to address the methods in the scaling function pointer arrays
                            // has to start from zero, so we have to map direction=-1 to di=0
                                         
    // Scale center
    length_type         scale_center_r;
    length_type         scale_center_z;                                               
    // Scale angle
    angle_type          scale_angle;
    Real                tan_scale_angle;
    
    // Array of scaling methods - we have to call different ones depending on the 
    // scale direction determined initially
    scaling_function_pt_type r1_function[2];
    scaling_function_pt_type z1_function[2];
    scaling_function_pt_type z2_function[2];
    
    // Important vectors and lengths used throughout the calculations
    // Local coordinate system in orthogonal scaling
    position_type       local_x, local_y, local_z;    
    position_type       ref_to_shell_vec;
    length_type         ref_to_shell_x, ref_to_shell_y, ref_to_shell_z;
    length_type         ref_to_shell_x2, ref_to_shell_y2, ref_to_shell_z2;
    length_type         scale_center_to_shell_x, scale_center_to_shell_y, scale_center_to_shell_z;
    
    
  public:       // CONSTRUCTOR
    
    CylinderScalingHelperTools( CylinderScalingFunctionsWrap<Ttraits_> *CSF_, 
                                position_type  testShell_reference_point_,
                                position_type  testShell_orientation_vector_,
                                position_type  testShell_dimensions_,
                                position_type  otherShell_position_t_,
                                position_type  otherShell_orientation_vector_,
                                position_type  otherShell_dimensions_,
                                direction_type direction_,
                                position_type  scale_center_info_, 
                                position_type  scale_angle_info_                 ):
     CSF(CSF_),
     testShell_reference_point(testShell_reference_point_),
     testShell_orientation_vector(testShell_orientation_vector_),
     testShell_dimensions(testShell_dimensions_),
     otherShell_position_t(otherShell_position_t_),
     otherShell_orientation_vector(otherShell_orientation_vector_),
     otherShell_dimensions(otherShell_dimensions_),
     direction(direction_),
     scale_center_info(scale_center_info_), 
     scale_angle_info(scale_angle_info_)
    {
         // Assign the scaling functions to the array for the two different directions
         // direction = -1 (di = 0)
         r1_function[0] = &CylinderScalingHelperTools<Ttraits_>::r_left;
         z1_function[0] = &CylinderScalingHelperTools<Ttraits_>::z_left;
         z2_function[0] = &CylinderScalingHelperTools<Ttraits_>::z_right;
         // direction = +1 (di = 1)
         r1_function[1] = &CylinderScalingHelperTools<Ttraits_>::r_right;
         z1_function[1] = &CylinderScalingHelperTools<Ttraits_>::z_right;
         z2_function[1] = &CylinderScalingHelperTools<Ttraits_>::z_left;
         
         // Set direction index for addressing the above arrays
         di = (direction == 1) ? 1 : 0;
         
         // Unpack scaling_info         
         scale_center_r  = scale_center_info[0];
         scale_center_z  = scale_center_info[1];         
         scale_angle     = scale_angle_info[0];
         tan_scale_angle = scale_angle_info[1];
         
         // Unpack shell dimensions
         // for testShell
         r  = testShell_dimensions[0];
         z1 = testShell_dimensions[1];
         z2 = testShell_dimensions[2];
         // for otherShell
         otherShell_radius = otherShell_dimensions[0];
         otherShell_hl     = otherShell_dimensions[1]; // half-length
         
         ref_to_shell_vec = subtract(otherShell_position_t, testShell_reference_point);
         ref_to_shell_z   = dot_product(ref_to_shell_vec, testShell_orientation_vector);
         
    };
    
    // DESTRUCTOR
    ~CylinderScalingHelperTools() {}; 

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
