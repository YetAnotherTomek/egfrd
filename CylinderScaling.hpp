#ifndef CYLINDERSCALING_HPP
#define CYLINDERSCALING_HPP

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>

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



// This is the structure that will hold the parameters for the h1- and r1-functions
// passed to the GSL rootfinder; we make it a template that will be parametrized
// from the CylinderScalingHelperTools class (with the same types as used there).
template <typename Ttraits_>
struct edge_hits_edge_eq_params
{ 
    typedef typename Ttraits_::length_type      length_type;
    typedef Real                                angle_type; // TODO could be improved
  
    angle_type  tan_scale_angle;
    length_type scale_center_to_shell_edge_x;
    length_type scale_center_to_shell_y;
    length_type scale_center_to_shell_z;
    length_type otherShell_radius_sq;
};

    
    

template <typename Ttraits_>
class CylinderScalingHelperTools
{
 
  typedef typename Ttraits_::length_type      length_type;  
  typedef typename Ttraits_::position_type    position_type;  
  typedef int                                 direction_type;
  typedef Real                                angle_type; // TODO could be improved
  
  enum collision_type { BARREL_HITS_FLAT, EDGE_HITS_EDGE, BARREL_HITS_EDGE, 
                        FLAT_HITS_BARREL, EDGE_HITS_BARREL, BARREL_HITS_BARREL };
                        
  enum collision_quadrant_type {
          LOWER_LEFT,   // otherShell overlaps with intersection of testShell in all four quadrants
          UPPER_LEFT,   // otherShell overlaps with intersection of testShell in the two quadrants closest to flat side of otherShell
          LOWER_RIGHT,  // otherShell overlaps with intersection of testShell in the two quadrants closest to projected barrel of otherShell
          UPPER_RIGHT,  // otherShell overlaps with intersection of testShell only in the quadrant closest to midpoint of otherShell          
          UNDETERMINED  // the default bogus value
       };
  
  // define scaling helper function reference type to later define an array
  // of scaling functions with direction_type index
  typedef length_type (CylinderScalingHelperTools<Ttraits_>::* scaling_function_pt_type)(length_type);
  // define reference to gsl function type
  typedef double (* gsl_function_pt_type)(double, void*);
  
  
  public:       // METHODS          
        
    /****************/
    /* MAIN METHOD  */    
    /****************/
    // Gives back the new cylinder dimensions.
    // Note that this only works properly if the whole class was correctly instantiated,
    // i.e. passing all necessary information via the constructor and storing it in
    // the required variables.
    inline position_type get_dr_dzright_dzleft_to_CylindricalShape()
    {
        // The master method that calls downstream helper methods in the right order
        // It returns the new testShell dimensions after scaling against otherShell,
        // after first determining whether the two shells are oriented in parallel 
        // or orthogonally
        
        length_type relative_orientation( dot_product(testShell_orientation_vector, otherShell_orientation_vector) );
        
        if( abs(feq(relative_orientation, 1.0, 1e-3 )) )
          return get_dr_dzright_dzleft_to_parallel_CylindricalShape();
        
        else if( abs(feq(relative_orientation, 0.0, 1e-3 )) )
          return get_dr_dzright_dzleft_to_orthogonal_CylindricalShape();
        
        else
          log_.warn("Shells seem neither parallel nor orthogonal in CylindricaltestShell scaling routine: relative_orientation = %.7f", relative_orientation);
          //throw unsupported("Shells are neither parallel nor orthogonal in CylindricaltestShell scaling routine."); // FIXME
        
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
    
    /*************************************/
    /* Collision of orthogonal cylinders */
    /*************************************/
    inline position_type get_dr_dzright_dzleft_to_orthogonal_CylindricalShape()
    {
        construct_local_coordinate_system();
        
        return get_dr_dzright_dzleft_for_specific_collision( determine_collision_type() );
    };
    
    /*********************************************************************************/
    /* Helper method to construct the local coordinate system for orthogonal scaling */
    /*********************************************************************************/
    inline void construct_local_coordinate_system()
    {
        // First define local_x depending on how the otherShell is oriented
        // with respect to the inter-shell-vector
        if(dot_product(ref_to_shell_vec, otherShell_orientation_vector) >= 0)
            local_x = otherShell_orientation_vector;
        else
            local_x = multiply(otherShell_orientation_vector, -1.0);
        
        // Now calculate local_z, taking into account on which side of 
        // testShell the otherShell is
        local_z = multiply(testShell_orientation_vector, direction);
        // local_y then is calculated as the cross product
        position_type cross( cross_product(local_x, local_z) );
        // However, we have some freedom of choice in the direction of local_y;
        // here we want it to point towards otherShell
        if(dot_product(ref_to_shell_vec, cross) >= 0)
            // cross already points towards the shell, we are fine
            local_y = cross;
        else
            // it points away from the shell, thus negate it
            local_y = multiply(cross, -1.0);

        // Transform the vector pointing from the ref. point towards the shell
        // into the local coordinate system
        ref_to_shell_x = dot_product(ref_to_shell_vec, local_x);
        ref_to_shell_y = dot_product(ref_to_shell_vec, local_y);
        ref_to_shell_z = dot_product(ref_to_shell_vec, local_z);
        assert((ref_to_shell_x >= 0.0) && (ref_to_shell_y >= 0.0) && (ref_to_shell_z >= 0.0));
              // by the definition of the local coordinate system
              
        // Also calculate the coordinates of the scale center
        // FIXME Attention: The scale center is not always on the z-axis!
        scale_center_to_shell_x = ref_to_shell_x;
        scale_center_to_shell_y = ref_to_shell_y;
        scale_center_to_shell_z = ref_to_shell_z - scale_center_z;
        // And the coordinates of the closest corner of a box surrounding the shape
        ref_to_shell_x2 = ref_to_shell_x - otherShell_hl;
        ref_to_shell_y2 = ref_to_shell_y - otherShell_radius;
        ref_to_shell_z2 = ref_to_shell_z - otherShell_radius; // not really needed, but defined here for completeness
        
    };
    
    /******************************************************************************/
    /* Helper method to determine how precisely the testShell hits the otherShell */
    /******************************************************************************/
    inline collision_type determine_collision_type()
    {
            
      collision_type collision_situation; // will be returned by this function      
      collision_quadrant_type quadrant = UNDETERMINED;
      
      if(ref_to_shell_x2 < 0.0 && ref_to_shell_y2 < 0.0)
      {
            // Quadrant 1
            quadrant = LOWER_LEFT;

            // Scale cylinder one until its height is equal to the z-distance minus the radius of the cyl. shell 2;
            // Then check whether any point of the axis of cyl. 2 projected on the top side of cyl. 1 is inside the radius of (scaled) cyl. 1.
            // If this is the case, we have the flat side of cyl. 1 hitting the barrel of cyl. 2. Else we have the edge hitting the barrel.
            length_type r1_touch( (scale_center_to_shell_z - otherShell_radius) * tan_scale_angle );

            if(r1_touch >= scale_center_to_shell_y)
                // At least a part of the projected axis is within r_touch => flat side of scaled cyl. hits the barrel of cyl. 2
                // Note that we know that at least a part of the shell is above the top side of the scaled cylinder (quadrant 1 condition)
                collision_situation = FLAT_HITS_BARREL;
            else
                collision_situation = EDGE_HITS_BARREL;
      }
      else if(ref_to_shell_x2 >= 0.0 && ref_to_shell_y2 < 0.0)
      {
            // Quadrant 2
            quadrant = UPPER_LEFT;

            length_type scale_center_to_shell_x_minus_half_length( scale_center_to_shell_x - otherShell_hl );
                        // just because this is reused several times below

            // The case scale_angle=0, i.e. radius remaining constant at scaling, has to be treated separately
            // because in this case the mathematics in the standard case misdetect the collision situation
            if(scale_angle == 0.0)
            {                
                if( std::sqrt( scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length
                               + scale_center_to_shell_y*scale_center_to_shell_y) < r )
                // The lowest point of the static cylinder is within the radius/flat side circle of the scaling cylinder ('testShell').
                // Therefore, when the height is scaled, the flat side must hit the barrel of the static cylinder at the lowpoint.
                    collision_situation = FLAT_HITS_BARREL;
                    
                else
                    collision_situation = EDGE_HITS_EDGE;
                    // This will also properly treat the case in which there is no collision because the projection of
                    // the static cylinder does not overlap with the flat side circle of the scaled cylinder
            }
            else
            {
                // case scale_angle >0

                // Two points on the edge of the static cylinder are of interest here:
                // - the "lowpoint", which is the point on the edge with the minimal z-distance to the scale center in the zy-plane
                // - the "critpoint" (critical point), which is the point on the edge with the shortest distance to the scale center                  
                length_type scale_center_to_flatend_x( scale_center_to_shell_x_minus_half_length - scale_center_r );
                length_type scale_center_to_lowpoint_x( std::sqrt( scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length
                                                                   + scale_center_to_shell_y*scale_center_to_shell_y ) - scale_center_r );
                                                        // TODO Is there no better way to include scale_center_r here?
                length_type scale_center_to_critpoint_z( scale_center_to_shell_z 
                                                         - std::sqrt(otherShell_radius*otherShell_radius - scale_center_to_shell_y*scale_center_to_shell_y) );
                length_type scale_center_to_lowpoint_z( scale_center_to_shell_z - otherShell_radius );

                // Strategy:
                // - When the scale angle is smaller than the angle between the scaled cylinder's axis and the line
                //   that links the critpoint then we have a BARREL_HITS_FLAT situation.
                // - When the scale angle is bigger than the angle between the scaled cylinder's axis and the line
                //   that links the lowpoint then we are in a FLAT_HITS_BARREL situation.                  
                // - Everything else results in EDGE_HITS_EDGE.

                // Calculate the angle between the scaled cylinder's axis and the line that links the scale center and lowpoint
                angle_type scale_center_to_shell_crit_angle_y( 0.0 );
                
                if(scale_center_to_critpoint_z == 0.0)
                    scale_center_to_shell_crit_angle_y = M_PI/2.0;                
                else
                {
                    scale_center_to_shell_crit_angle_y = std::atan(scale_center_to_flatend_x/ scale_center_to_critpoint_z);
                    
                    if(scale_center_to_critpoint_z < 0.0)
                        scale_center_to_shell_crit_angle_y += M_PI;
                }

                // Calculate the angle between the scaled cylinder's axis and the line that links the scale center and critpoint
                angle_type scale_center_to_shell_low_angle_y( 0.0 );
                
                if(scale_center_to_lowpoint_z == 0.0)
                    scale_center_to_shell_low_angle_y = M_PI/2.0;
                else
                {
                    scale_center_to_shell_low_angle_y  = std::atan(scale_center_to_lowpoint_x/ scale_center_to_lowpoint_z );
                    if(scale_center_to_lowpoint_z < 0.0)
                        scale_center_to_shell_low_angle_y += M_PI;
                }

                // Compare the angles to determine the collision situation
                if(scale_angle <= scale_center_to_shell_crit_angle_y)
                    collision_situation = BARREL_HITS_FLAT;
                else if(scale_angle >= scale_center_to_shell_low_angle_y)
                    collision_situation = FLAT_HITS_BARREL;
                else
                {
                    assert(scale_center_to_shell_crit_angle_y < scale_angle
                              && scale_angle < scale_center_to_shell_low_angle_y);
                    collision_situation = EDGE_HITS_EDGE;
                    // r1_min and h1_min are later used in the EDGE_HITS_EDGE rootfinder-based collision routine
                    // TODO can we not outsource this from here?
                    this->r1_min = (scale_center_to_shell_x - otherShell_hl)*(1.0+TOLERANCE);
                    this->h1_min = this->r1_min / tan_scale_angle;
                }
            }
      }
      else if(ref_to_shell_x2 < 0.0 && ref_to_shell_y2 >= 0.0)
      {
            // Quadrant 3
            quadrant = LOWER_RIGHT;

            // The case scale_angle=0, i.e. radius remaining constant at scaling, has to be treated separately
            // because in this case the mathematics in the standard case misdetect the collision situation
            if(scale_angle == 0.0)
            {
                if(scale_center_to_shell_y < r)
                // The lowest point of the static cylinder is within the radius/flat side circle of the scaling cylinder ('testShell').
                // Therefore, when the height is scaled, the flat side must hit the barrel of the static cylinder at the lowpoint.
                    collision_situation = FLAT_HITS_BARREL;

                else if(scale_center_to_shell_y - otherShell_radius <= r)
                    collision_situation = EDGE_HITS_BARREL;

                else
                    // In this case the cylinders do not hit; this is treated properly by the BARREL_HITS_BARREL routine
                    collision_situation = BARREL_HITS_BARREL;
            }
            else
            {
                // Two points on the edge of the static cylinder are of interest here:
                // - the "lowpoint", which is the point on the edge with the minimal z-distance to the scale center in the zy-plane
                // - the "critpoint" (critical point), which is the point on the edge with the shortest distance to the scale center   
                scale_center_to_shell_y -= scale_center_r;  // TODO Is there no better way to include scale_center_r here?
                                                            // This is dangerous, permanently changing this.scale_center_to_shell_y!
                
                length_type scale_center_to_critpoint_y( scale_center_to_shell_y - otherShell_radius );
                length_type scale_center_to_lowpoint_z(  scale_center_to_shell_z - otherShell_radius );

                // Calculate the angle between the scaled cylinder's axis and the line that links the scale center and critpoint
                angle_type scale_center_to_shell_crit_angle_x( 0.0 );
                if(scale_center_to_shell_z == 0.0)
                    scale_center_to_shell_crit_angle_x = M_PI/2.0;
                else
                {
                    scale_center_to_shell_crit_angle_x = std::atan(scale_center_to_critpoint_y / scale_center_to_shell_z);
                    if(scale_center_to_shell_z < 0.0)
                        scale_center_to_shell_crit_angle_x += M_PI;
                }
                // Calculate the angle between the scaled cylinder's axis and the line that links the scale center and lowpoint
                angle_type scale_center_to_shell_low_angle_x( 0.0 );
                if(scale_center_to_lowpoint_z == 0.0)
                    scale_center_to_shell_low_angle_x = M_PI/2.0;
                else
                {
                    scale_center_to_shell_low_angle_x  = std::atan(scale_center_to_shell_y / scale_center_to_lowpoint_z);
                    if(scale_center_to_lowpoint_z < 0.0)
                        scale_center_to_shell_low_angle_x += M_PI;
                }
                // Determine which collision will happen:
                // If the scale angle is bigger than the angle between the scaled cylinder's axis and the lowpoint
                // then we hit the barrel of the static cylinder when scaling.
                // If the scale angle is smaller than the angle between the scaled cylinder's axis and the critpoint
                // the barrel of the scaled cylinder will hit the barrel of the static cylinder.
                // Otherwise we have the edge of the scaled cylinder hitting the barrel of the static cylinder.
                if(scale_angle <= scale_center_to_shell_crit_angle_x)
                    collision_situation = BARREL_HITS_BARREL;
                else if(scale_center_to_shell_low_angle_x <= scale_angle)
                    collision_situation = FLAT_HITS_BARREL;
                else
                {
                    assert(scale_center_to_shell_crit_angle_x < scale_angle
                              && scale_angle < scale_center_to_shell_low_angle_x);
                    collision_situation = EDGE_HITS_BARREL;
                }
            }
      }
      else
      {
            // Quadrant 4
            quadrant = UPPER_RIGHT;

            assert(ref_to_shell_x2 >= 0.0 && ref_to_shell_y2 >= 0.0); // this just must be true!

            length_type scale_center_to_shell_x_minus_half_length( scale_center_to_shell_x - otherShell_hl );
            length_type scale_center_to_shell_x_minus_half_length_sq( scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length );
            length_type scale_center_to_shell_y_minus_radius( scale_center_to_shell_y - otherShell_radius );

            // Two points on the edge of the static cylinder are of interest here:
            // - the "lowpoint", which is the point on the edge with the minimal z-distance to the scale center in the zy-plane
            // - the "critpoint" (critical point), which is the point on the edge with the minimal y-distance to
            //   the scale center in the zy-plane
            length_type scale_center_to_critpoint_r( std::sqrt(scale_center_to_shell_y_minus_radius*scale_center_to_shell_y_minus_radius
                                                               + scale_center_to_shell_x_minus_half_length_sq) - scale_center_r );  // a_r
            length_type scale_center_to_lowpoint_r( std::sqrt(scale_center_to_shell_y*scale_center_to_shell_y 
                                                              + scale_center_to_shell_x_minus_half_length_sq) - scale_center_r );  // b_r
            length_type scale_center_to_lowpoint_z( scale_center_to_shell_z - otherShell_radius );  // b_z

            // Calculate the angle between the z-axis of the scaled cylinder and the line that
            // connects this axis with the "critpoint" on the static cylinder's edge
            angle_type scale_center_to_shell_crit_angle_xy( 0.0 );
            if(scale_center_to_shell_z == 0.0)
                scale_center_to_shell_crit_angle_xy = M_PI/2.0;
            else
            {
                scale_center_to_shell_crit_angle_xy = std::atan(scale_center_to_critpoint_r / scale_center_to_shell_z);
                if(scale_center_to_shell_z < 0.0)
                    scale_center_to_shell_crit_angle_xy += M_PI;
            }

            // Calculate the angle between the z-axis of the scaled cylinder and the line that
            // connects this axis with the "lowpoint" on the static cylinder's edge
            angle_type scale_center_to_shell_low_angle_xy( 0.0 );
            if(scale_center_to_lowpoint_z == 0.0)
                scale_center_to_shell_low_angle_xy = M_PI/2.0;
            else
            {
                scale_center_to_shell_low_angle_xy  = std::atan(scale_center_to_lowpoint_r / scale_center_to_lowpoint_z );
                if(scale_center_to_lowpoint_z < 0.0)
                    scale_center_to_shell_low_angle_xy += M_PI;
            }

            // First treat the special case: when the scale center is further away from the z-axis than the "critpoint" (in the xy-plane).
            // This also treats the case in which scale_angle = 0 and the scaled cylinder is directly below the static shell.
            if(scale_center_to_critpoint_r <= 0.0)
            {
                collision_situation = EDGE_HITS_EDGE;

                if(scale_center_to_lowpoint_r <= 0.0)
                    collision_situation = FLAT_HITS_BARREL;
            }
            // Now the standard case:
            else
            {
                if(scale_angle <= scale_center_to_shell_crit_angle_xy)
                    collision_situation = BARREL_HITS_EDGE;

                else if(scale_center_to_shell_low_angle_xy <= scale_angle)
                    collision_situation = FLAT_HITS_BARREL;

                else
                {
                    assert(scale_center_to_shell_crit_angle_xy < scale_angle 
                           && scale_angle < scale_center_to_shell_low_angle_xy 
                           && scale_angle > 0.0);
                           
                    collision_situation = EDGE_HITS_EDGE;
                    // r1_min and h1_min are later used in the EDGE_HITS_EDGE rootfinder-based collision routine 
                    // TODO can we not outsource this from here?
                    this->r1_min = std::sqrt(scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length +
                                             scale_center_to_shell_y_minus_radius*scale_center_to_shell_y_minus_radius) * (1.0+TOLERANCE);
                    this->h1_min = this->r1_min / tan_scale_angle;
                }
            }
      }
      
      return collision_situation;
      
    };
    
    /*********************************************************************************************/
    /* Helper method that calculates the new lengths for the specific collision determined above */
    /*********************************************************************************************/
    // This function partly uses a rootfinder method to determine the new test shell dimensions
    // (in particular in the case in which the two orthogonal cylinder's edges hit each other).
    // The associated helper scaling functions and a common parameter holder are defined first here:   
    
    // (a) Helper function to find h1 as a function of r1 in edge_hits_edge scaling
    static length_type edge_hits_edge_h1_eq(length_type x, void* params)
    {
        log_.info("     FROM h1 FUNCTION:");
        // Cast parameters into the right form
        struct edge_hits_edge_eq_params<Ttraits_>* p = (struct edge_hits_edge_eq_params<Ttraits_>*) params;
        // Unpack parameters
        log_.info("     Unpacking...");
        Real         tan_scale_angle              = (p->tan_scale_angle);
        length_type  scale_center_to_shell_edge_x = (p->scale_center_to_shell_edge_x);
        length_type  scale_center_to_shell_y      = (p->scale_center_to_shell_y);
        length_type  scale_center_to_shell_z      = (p->scale_center_to_shell_z);
        length_type  otherShell_radius_sq         = (p->otherShell_radius_sq);
                
        log_.info("     tan_scale_angle              = %.20f", tan_scale_angle);
        log_.info("     scale_center_to_shell_edge_x = %.20f", scale_center_to_shell_edge_x);
        log_.info("     scale_center_to_shell_y      = %.20f", scale_center_to_shell_y);
        log_.info("     scale_center_to_shell_z      = %.20f", scale_center_to_shell_z);
        log_.info("     otherShell_radius_sq         = %.20f", otherShell_radius_sq);
      
        // We will take the square root of the following below
        Real sqrt_arg( (x*tan_scale_angle)*(x*tan_scale_angle) - scale_center_to_shell_edge_x*scale_center_to_shell_edge_x );

        if(sqrt_arg < 0.0 && std::abs(sqrt_arg) <= TOLERANCE*scale_center_to_shell_edge_x*scale_center_to_shell_edge_x)
        {
                sqrt_arg = 0.0;     // This safety check is to prevent math domain errors
                                    // in case sqrt_arg is close to zero and taking the
                                    // difference results in very small negative numbers
                log_.warn("Orthogonal cylinder scaling, EDGE_HITS_EDGE case: Setting small negative sqrt argument to zero within TOLERANCE.");
        }
                          
        return (scale_center_to_shell_z - x)*(scale_center_to_shell_z - x) - otherShell_radius_sq 
                + (scale_center_to_shell_y - std::sqrt( sqrt_arg ))*(scale_center_to_shell_y - std::sqrt( sqrt_arg ));
    };
    
    // (b) Helper function to find r1 as a function of h1 in edge_hits_edge scaling
    // Very similar, but in a crucial part different from edge_hits_edge_h1_eq() defined above
    static length_type edge_hits_edge_r1_eq(length_type x, void* params)
    {   
        log_.info("     FROM r1 FUNCTION:");
        // Cast parameters into the right form
        struct edge_hits_edge_eq_params<Ttraits_>* p = (struct edge_hits_edge_eq_params<Ttraits_>*) params;        
        // Unpack parameters
        log_.info("     Unpacking...");
        Real         tan_scale_angle              = (p->tan_scale_angle);
        length_type  scale_center_to_shell_edge_x = (p->scale_center_to_shell_edge_x);
        length_type  scale_center_to_shell_y      = (p->scale_center_to_shell_y);
        length_type  scale_center_to_shell_z      = (p->scale_center_to_shell_z);
        length_type  otherShell_radius_sq         = (p->otherShell_radius_sq);
        
        assert(tan_scale_angle != 0.0);
                
        log_.info("     tan_scale_angle              = %.20f", tan_scale_angle);
        log_.info("     scale_center_to_shell_edge_x = %.20f", scale_center_to_shell_edge_x);
        log_.info("     scale_center_to_shell_y      = %.20f", scale_center_to_shell_y);
        log_.info("     scale_center_to_shell_z      = %.20f", scale_center_to_shell_z);
        log_.info("     otherShell_radius_sq         = %.20f", otherShell_radius_sq);
      
        // We will take the square root of the following below
        Real sqrt_arg( x*x - scale_center_to_shell_edge_x*scale_center_to_shell_edge_x );
              
        if(sqrt_arg < 0.0 && std::abs(sqrt_arg) <= TOLERANCE*scale_center_to_shell_edge_x*scale_center_to_shell_edge_x)
        {
              sqrt_arg = 0.0;     // This safety check is to prevent math domain errors
                                  // in case sqrt_arg is close to zero and taking the
                                  // difference results in very small negative numbers
              log_.warn("Orthogonal cylinder scaling, EDGE_HITS_EDGE case: Setting small negative sqrt argument to zero within TOLERANCE.");
        }

        return (scale_center_to_shell_z - x/tan_scale_angle)*(scale_center_to_shell_z - x/tan_scale_angle) - otherShell_radius_sq
                + (scale_center_to_shell_y - std::sqrt( sqrt_arg ))*(scale_center_to_shell_y - std::sqrt( sqrt_arg ));
    };
               
    // (c) Now the function that handles the scaling for each collision situation specifically
    inline position_type get_dr_dzright_dzleft_for_specific_collision(collision_type collision_situation)
    {
      
      // The return lengths
      length_type r_new, z1_new, z2_new; // TODO what happens with z2_new?
      
      // Some commonly used values
      length_type otherShell_radius_sq( otherShell_radius*otherShell_radius );
      
      log_.info("C++: Handling specific collision situation...");  // TESTING
      
      switch(collision_situation)
      {
            
          case BARREL_HITS_FLAT:
          { 
              // the scaling cylinder ('testShell') hits the flat side of 'otherShell' with its barrel side
              r_new  = std::min(r, (ref_to_shell_x - otherShell_hl));
              z1_new = std::min(z1, (this->*z1_function[di])(r_new));
          }
          break;
          
          case EDGE_HITS_EDGE:
          {      
              log_.info("C++: Entering EDGE_HITS_EDGE case...");  // TESTING
              // the scaling cylinder ('testShell') hits the edge of 'otherShell' with its edge
              // TODO we have a solution but it can only be found with a root finder -> slow
              
              // We need the distances to the edge of otherShell in the following:
              length_type scale_center_to_shell_edge_x( scale_center_to_shell_x - otherShell_hl );
              length_type scale_center_to_shell_edge_y( scale_center_to_shell_y - otherShell_radius );
              length_type scale_center_to_shell_edge_z( scale_center_to_shell_z - otherShell_radius );
              // Also pre-calculate the squares
              length_type scale_center_to_shell_edge_x_sq( scale_center_to_shell_edge_x*scale_center_to_shell_edge_x );
              length_type scale_center_to_shell_edge_y_sq( scale_center_to_shell_edge_y*scale_center_to_shell_edge_y );
              length_type scale_center_to_shell_edge_z_sq( scale_center_to_shell_edge_z*scale_center_to_shell_edge_z );
              // Also used below
              length_type scale_center_to_shell_y_sq( scale_center_to_shell_y*scale_center_to_shell_y );


              if( scale_angle == 0.0 )
              {
                  // scale_angle = 0 means that the radius r will stay constant at scaling
                  // Therefore, if the projected edge of the static shell lies outside of the 
                  // circle defined by r there will be no collision. We check for that:
                  if(r < std::sqrt( scale_center_to_shell_edge_x_sq + scale_center_to_shell_edge_y_sq))
                  {
                      // The height and radius of the scaled cylinder remain unaltered
                      r_new = r;
                      z1_new = z1;        // Note that z1_function would give infinity,
                                          // but we would take the min(z1_function(r_new), z1)
                  }
                  else
                  {
                      // There is a collision possible when z1 is scaled                    
                      length_type y1_collide( std::sqrt(r*r - scale_center_to_shell_edge_x_sq) );
                      length_type y2_collide( scale_center_to_shell_y - y1_collide );
                      length_type h_collide(  std::sqrt(otherShell_radius_sq - y2_collide*y2_collide) );

                      length_type h_touch( scale_center_to_shell_z - h_collide );
                      z1_new = std::min(z1, h_touch);
                      r_new  = std::min(r,  (this->*r1_function[di])(z1_new));
                  }
              }
              else
              {
                  // TESTING Code snippet to check whether rootfinder run is necessary
                  // We only want to run it if the radius increase at intersection of the cylinder edges is "significant".
                  // This is sub-optimal but it avoids rootfinder errors caused by the fact that the radius search interval
                  // becomes to small, which is equivalent to the potential r-increase becoming small.
                  //
                  // We define the increase as "significant" by comparing the potential radius increase with the radial (xy-) distance
                  // from the scale center to the static shell; if increase / distance >= 10% the gain is "significant".
                  // adial_distance = math.sqrt( scale_center_to_shell_x**2 + scale_center_to_shell_y**2 )
                  // potential_r_increase = math.sqrt( scale_center_to_shell_y**2 + scale_center_to_shell_edge_x**2 ) - scale_center_to_shell_edge_x  
                  // assert radial_distance >= potential_r_increase and potential_r_increase >= 0
                  /*
                      bool run_rootfinder( true );
                      if( potential_r_increase / radial_distance >= 0.10 )
                          run_rootfinder = true;
                      else
                      {
                          run_rootfinder = false;
              
                          length_type r_touch( scale_center_to_shell_edge_y );
                          if(scale_center_to_shell_x < scale_center_to_shell_y)
                                  r_touch = scale_center_to_shell_edge_y;                    
                          
                          r_new  = std::min(r, r_touch);
                          z1_new = std::min(z1, (this->*z1_function[di])(r_new) );
                      }

                      //if(run_rootfinder) ...
                  */

                  // In this case there is no analytical solution, we have to use a rootfinder.
                  // To this purpuse we make use of the GSL Brent rootfinder.

                  // Construct a rootfinder instance based on the function defined above                  
                  log_.info("C++: Defining rootfinder and parameter structure");  // TESTING
                  const gsl_root_fsolver_type* solver;
                  gsl_root_fsolver* solver_ref;
                  int rf_status;
                  int iter = 0, max_iter = 100;
                  double lo, hi;
                  double root;
                  // The function to be solved has to be casted to the right (GSL) format
                  gsl_function F;
                  // The parameters are required to be passed by a grouping structure
                  edge_hits_edge_eq_params<Ttraits_> p = { tan_scale_angle, 
                                                           scale_center_to_shell_edge_x, 
                                                           scale_center_to_shell_y, 
                                                           scale_center_to_shell_z,
                                                           otherShell_radius_sq
                                                         };
                  log_.info("     PARAMETERS STRUCTURE:"); // TESTING
                  log_.info("     scale_angle                  = %.20f", (double)scale_angle);
                  log_.info("     tan_scale_angle              = %.20f", (double)tan_scale_angle);
                  log_.info("     scale_center_to_shell_edge_x = %.20f", (double)scale_center_to_shell_edge_x);
                  log_.info("     scale_center_to_shell_y      = %.20f", (double)scale_center_to_shell_y);
                  log_.info("     scale_center_to_shell_z      = %.20f", (double)scale_center_to_shell_z);
                  log_.info("     otherShell_radius_sq         = %.20f", (double)otherShell_radius_sq);
                    
                  
                  if(scale_angle <= M_PI/4.0)
                  {                    
                   
                    log_.info("C++: Entering rootfinding, scale_angle <= M_PI/4.0 (finding h1)");  // TESTING
                                                                            
                    // Pack function and parameters into the format required by GSL
                    F.function = reinterpret_cast<typeof(F.function)>( &CylinderScalingHelperTools<Ttraits_>::edge_hits_edge_h1_eq );
                    F.params = &p;
                    
                    // TESTING call
                    log_.info("C++: Rootfinder function test call:");
                    //(*F.function)(1000000.0, &p);
                    //(*F.function)(1000000.0, F.params);
                    
                    // Set the iteration bounds
                    length_type h1_interval_min( (this->*z1_function[di])(scale_center_to_shell_edge_x) - scale_center_z );
                    length_type h1_interval_max( (this->*z1_function[di])( std::sqrt( scale_center_to_shell_y_sq + scale_center_to_shell_edge_x_sq ) )
                                                  - scale_center_z );
                                                  
                    assert(h1_interval_min >= 0.0);
                    assert(h1_interval_max >= h1_interval_min);
                                                        
                    // Create a rootfinder instance with pointer
                    solver     = gsl_root_fsolver_brent;
                    solver_ref = gsl_root_fsolver_alloc(solver);
                    // Initialize
                    log_.info("C++: Initializing rootfinder.");  // TESTING                        
                    gsl_root_fsolver_set(solver_ref, &F, h1_interval_min, h1_interval_max);
                    
                    log_.info("C++: Starting rootfinder iteration."); // TESTING                    
                    root = z1;
                    iter = 0;
                    do
                    {
                        iter++;
                        rf_status = gsl_root_fsolver_iterate(solver_ref);
                        root = gsl_root_fsolver_root(solver_ref);
                        lo   = gsl_root_fsolver_x_lower(solver_ref);
                        hi   = gsl_root_fsolver_x_upper(solver_ref);
                        rf_status = gsl_root_test_interval(lo, hi, 0, 0.001); // arg. 3 is absolute error, arg. 4 relative error of the found root
                  
                        if (rf_status == GSL_SUCCESS)
                              log_.info ("     Converged:");
                  
                        log_.info("     Iteration %4d of %4d: bounds=[%.7f, %.7f], value=%.7f, error=%.7f",
                                                  iter, max_iter, lo, hi, root, hi - lo);
                    }
                    while (rf_status == GSL_CONTINUE && iter < max_iter);                                    
                    
                    // Finally, assemble the solution; z1_new and r_new will be passed back by this function further below
                    length_type h_touch( scale_center_z + root );
                    z1_new = std::min(z1, h_touch);
                    r_new  = std::min(r,  (this->*r1_function[di])(z1_new));
                    
                    // Cleanup
                    gsl_root_fsolver_free(solver_ref);


                  // TODO TODO TODO - Maybe we don't need that iterative interval adaptation magic any more with the GSL rootfinder
                  /* 
                      // To ensure that we always search on an interval within which h1_eq() contains only one zero root
                      // we enlarge the rootfinder interval, starting from safe bounds, towards the bounds beyond which we
                      // know there are unwanted solutions. We stop the enlargement when function h1_eq() certainly 
                      // changes sign within the interval.
                      // We use a multiplicative factor tau to enlarge the domain. tau is set such that 
                      // (1-tau)*h1_interval_max - (1+tau)*h1_interval_min >= 0.5*(h1_interval_max-h1_interval_min)
                      // in order to ensure that the interval is nonzero and positive in the first iteration.

                      // Set some repeating and default values
                      // r1_interval min, r1_interval_max are the outer interval bounds for which there is only one root
                      h1_interval_min = (this->*z1_function[di])(scale_center_to_shell_edge_x) - scale_center_z
                      h1_interval_max = (this->*z1_function[di])( math.sqrt( scale_center_to_shell_y*scale_center_to_shell_y + \
                                                        scale_center_to_shell_edge_x*scale_center_to_shell_edge_x ) ) - scale_center_z
                      h1_interval_start = h1_interval_min
                      h1_interval_end   = h1_interval_max
                      // Set the enlargement factor; we want it to be at least as small as TOLERANCE
                      tau = min(TOLERANCE, 0.5*(h1_interval_max-h1_interval_min)/(h1_interval_max+h1_interval_min))
                      assert tau>0.0 and tau<1.0

                      // Start the iteration
                      n=1
                      nmax=100
                        // Since tau is very small, the iteration hardly ever should get to nmax; but to be safe...
                      h1_eq_product = 1                    
                      while h1_eq_product > 0.0:

                          h1_interval_start = (1.0+tau**n) * h1_interval_min
                          h1_interval_end   = (1.0-tau**n) * h1_interval_max
                          h1_eq_product     = h1_eq(h1_interval_start) * h1_eq(h1_interval_end)

                          n = n+1

                          if n>nmax:
                              raise testShellError('get_dr_dzright_dzleft_to_CylindricalShape: Could not find suitable rootfinder boundaries in EDGE_HITS_EDGE cylinder scaling. nmax=%s' % str(nmax) )

                      assert(h1_interval_start >= 0)
                      assert(h1_interval_end   >= h1_interval_end)

                      // TODO TESTING REMOVE THIS WHEN DONE
                      // log.debug( "***** NEW ROOTFINDER ITERATION *****" )
                      // log.debug( "  Dy-r2 = %s" % (scale_center_to_shell_y-otherShell_radius) )
                      // log.debug( "  Dx-h2 = %s" % scale_center_to_shell_edge_x )
                      // log.debug( "  h1_interval_start = %s" % h1_interval_start )
                      // log.debug( "  h1_interval_end = %s"   % h1_interval_end )
                      // log.debug( "  h1(i_start) = %s"       % h1_eq(h1_interval_start) )
                      // log.debug( "  h1(i_end) = %s"         % h1_eq(h1_interval_end) )

                      // Finally, start the rootfinding
                      h_touch = scale_center_z + findroot(h1_eq, h1_interval_start, h1_interval_end)
                      z1_new = std::min(z1, h_touch)
                      r_new  = std::min(r,  (this->*r1_function[di])(z1_new))
                      
                  */ // TODO TODO TODO
                  
                  }
                  else // if scale_angle > M_PI/4.0
                  {
                    log_.info("C++: Entering rootfinding, scale_angle > M_PI/4.0 (finding r1)");  // TESTING                                        
                                                        
                    // Pack function and parameters into the format required by GSL
                    F.function = reinterpret_cast<typeof(F.function)>( &CylinderScalingHelperTools<Ttraits_>::edge_hits_edge_r1_eq );
                    F.params = &p;
                    
                    // TESTING call
                    log_.info("C++: Rootfinder function test call:");
                    (*F.function)(1000000.0, &p);
                    (*F.function)(1000000.0, F.params);
                    
                    // Set the iteration bounds
                    length_type r1_interval_min( scale_center_to_shell_edge_x );
                    length_type r1_interval_max( std::sqrt( scale_center_to_shell_y_sq + scale_center_to_shell_edge_x_sq ) );
                    
                    assert(r1_interval_min >= 0.0);
                    assert(r1_interval_max >= r1_interval_min);
                                                                          
                    // Create a rootfinder instance with pointer
                    solver     = gsl_root_fsolver_brent;
                    solver_ref = gsl_root_fsolver_alloc(solver);
                    // Initialize
                    log_.info("C++: Initializing rootfinder.");  // TESTING
                    gsl_root_fsolver_set(solver_ref, &F, r1_interval_min, r1_interval_max);
                    
                    log_.info("C++: Starting rootfinder iteration."); // TESTING
                    root = r;
                    iter = 0;
                    do
                    {
                        iter++;
                        rf_status = gsl_root_fsolver_iterate(solver_ref);
                        root = gsl_root_fsolver_root(solver_ref);
                        lo   = gsl_root_fsolver_x_lower(solver_ref);
                        hi   = gsl_root_fsolver_x_upper(solver_ref);
                        rf_status = gsl_root_test_interval(lo, hi, 0, 0.001); // arg. 3 is absolute error, arg. 4 relative error of the found root
                  
                        if (rf_status == GSL_SUCCESS)
                              log_.info ("     Converged:");
                  
                        log_.info("     Iteration %4d of %4d: bounds=[%.7f, %.7f], value=%.7f, error=%.7f",
                                                  iter, max_iter, lo, hi, root, hi - lo);
                    }
                    while (rf_status == GSL_CONTINUE && iter < max_iter);                
                    
                    // Finally, assemble the solution; z1_new and r_new will be passed back by this function further below                    
                    r_new  = std::min(r, root);
                    z1_new = std::min(z1, (this->*z1_function[di])(r_new));
                    log_.info("C++: Result after rootfinding: r_new=%.7f, z1_new=%.7f", r_new, z1_new);
                    
                    // Cleanup
                    gsl_root_fsolver_free(solver_ref);
                    
                    
                    /* // TODO TODO TODO
                    
                      // This uses the same rootfinding equation and interval as above with r1=h1/tan_scale_angle
                      // instead of h1; notice that scale_angle > 0 here, so we can divide by it
                      def r1_eq(x):

                          // We will take the square root of the following below
                          sqrt_arg = x*x - scale_center_to_shell_edge_x*scale_center_to_shell_edge_x

                          if sqrt_arg < 0.0 and abs(sqrt_arg) <= TOLERANCE*scale_center_to_shell_edge_x*scale_center_to_shell_edge_x:

                              sqrt_arg = 0.0      // This safety check is to prevent math domain errors
                                                  // in case sqrt_arg is close to zero and taking the
                                                  // difference results in very small negative numbers
                              log.warn('Orthogonal cylinder scaling, EDGE_HITS_EDGE case: Setting negative sqrt argument to zero within TOLERANCE.')

                          eq_value = (scale_center_to_shell_z - x/tan_scale_angle)*(scale_center_to_shell_z - x/tan_scale_angle) - otherShell_radius_sq +\
                                        (scale_center_to_shell_y - math.sqrt( sqrt_arg ) )*(scale_center_to_shell_y - math.sqrt( sqrt_arg ) )

                          // TODO TESTING REMOVE THIS WHEN DONE
                          // log.debug( "***** ROOTFINDER CALL: Calling r1_eq() with: *****" )
                          // log.debug( "  r1 = %s, value = %s" % (x, eq_value ) )
                          // log.debug( "  h1 = %s" % ((this->*z1_function[di])(x)) )

                          return eq_value

                      // To ensure that we always search on an interval within which r1_eq() contains only one zero root
                      // we enlarge the rootfinder interval, starting from safe bounds, towards the bounds beyond which we
                      // know there are unwanted solutions. We stop the enlargement when function r1_eq() certainly 
                      // changes sign within the interval.
                      // We use a multiplicative factor tau to enlarge the domain. tau is set such that 
                      // (1-tau)*r1_interval_max - (1+tau)*r1_interval_min >= 0.5*(r1_interval_max-r1_interval_min)
                      // in order to ensure that the interval is nonzero and positive in the first iteration.

                      // Set some repeating and default values
                      // r1_interval min, r1_interval_max are the outer interval bounds for which there is only one root
                      r1_interval_min = scale_center_to_shell_edge_x
                      r1_interval_max = math.sqrt( scale_center_to_shell_y*scale_center_to_shell_y + \
                                                      scale_center_to_shell_edge_x*scale_center_to_shell_edge_x )
                      r1_interval_start = r1_interval_min
                      r1_interval_end   = r1_interval_max
                      // Set the enlargement factor; we want it to be at least as small as TOLERANCE
                      tau = min(TOLERANCE, 0.5*(r1_interval_max-r1_interval_min)/(r1_interval_max+r1_interval_min))
                      assert tau>0.0 and tau<1.0

                      // Start the iteration
                      n=1
                      nmax=100
                        // Since tau is very small, the iteration hardly ever should get to nmax; but to be safe...
                      r1_eq_product = 1
                      while r1_eq_product > 0.0:

                          r1_interval_start = (1.0+tau**n) * r1_interval_min
                          r1_interval_end   = (1.0-tau**n) * r1_interval_max
                          r1_eq_product     = r1_eq(r1_interval_start) * r1_eq(r1_interval_end)
                          
                          n = n+1

                          if n>nmax:
                              raise testShellError('get_dr_dzright_dzleft_to_CylindricalShape: Could not find suitable rootfinder boundaries in EDGE_HITS_EDGE cylinder scaling. nmax=%s' % str(nmax) )

                      // Non-adaptive version; TODO Remove this after TESTING
                      // r1_interval_start = scale_center_to_shell_edge_x
                      // r1_interval_end   = math.sqrt( scale_center_to_shell_y**2 + scale_center_to_shell_edge_x**2 )
                      assert(r1_interval_start >= 0), 'r1_interval_start = %s, n_iterations=%s' % (r1_interval_start, n)
                      assert(r1_interval_end   >= r1_interval_start), 'r1_interval_start = %s, r1_interval_end=%s, n_iterations=%s' % (r1_interval_start, r1_interval_end, n)

                      // TODO TESTING REMOVE THIS WHEN DONE
                      // print "***** NEW ROOTFINDER ITERATION *****"
                      // print "  Dy-r2 = %s" % (scale_center_to_shell_y-otherShell_radius)
                      // print "  Dx-h2 = %s" % scale_center_to_shell_edge_x
                      // print "  r1_interval_start = %s" % r1_interval_start
                      // print "  r1_interval_end = %s"   % r1_interval_end
                      // print "  r1(i_start) = %s"       % r1_eq(scale_center_to_shell_edge_x)
                      // print "  r1(i_end) = %s"         % r1_eq(math.sqrt( scale_center_to_shell_y**2 + scale_center_to_shell_edge_x**2))

                      // Finally, start the rootfinding
                      r_touch = findroot(r1_eq, r1_interval_start, r1_interval_end)
                                
                      r_new  = min(r, r_touch)
                      z1_new = min(z1, (this->*z1_function[di])(r_new))
                  */ // TODO TODO TODO                  
                  }                 
                  
              }; // if scale_angle <> 0.0                                          
              
          }
          break;
          
          case BARREL_HITS_EDGE:
          {   
              // the scaling cylinder ('testShell') hits the edge of 'otherShell' with its barrel side
              r_new = std::min(r, std::sqrt( (ref_to_shell_x-otherShell_hl)*(ref_to_shell_x-otherShell_hl) \
                                              + (ref_to_shell_y - otherShell_radius)*(ref_to_shell_y - otherShell_radius) ));
              z1_new = std::min(z1, (this->*z1_function[di])(r_new));                  
          }
          break;

          case FLAT_HITS_BARREL:
          {   
              // the scaling cylinder ('testShell') hits the barrel of 'otherShell' with its top flat side
              z1_new = std::min(z1, (ref_to_shell_z - otherShell_radius));
              r_new  = std::min(r,  (this->*r1_function[di])(z1_new));
          }     
          break;

          case EDGE_HITS_BARREL:
          {
              // the scaling cylinder ('testShell') hits the barrel of 'otherShell' with its edge
              
              // If scale_angle == 0, the scale center is offset from the reference point,
              // located at the barrel of the scaled cylinder. In order for the below calculation
              // to succeed we have to take this into account. In general scale_center_r == 0
              // so that this step only affects the scale_angle == 0 cases:
              scale_center_to_shell_y -= scale_center_r;

              // Pre-calculate some lengths...
              length_type ss_sq( scale_center_to_shell_z*scale_center_to_shell_z + scale_center_to_shell_y*scale_center_to_shell_y );
              length_type scale_center_to_shell( std::sqrt(ss_sq) );
              // ...and angles used below
              angle_type shell_angle_yz( std::atan(scale_center_to_shell_y/scale_center_to_shell_z) );
              angle_type angle_diff( std::abs(shell_angle_yz - scale_angle) );
              Real sin_angle_diff( std::sin(angle_diff) );
              Real cos_angle_diff( std::cos(angle_diff) ); // TODO fix typing
              
              assert( scale_center_to_shell >= otherShell_radius );// should never fail

              angle_type ss_angle( M_PI - std::asin(std::sin(angle_diff)*scale_center_to_shell/otherShell_radius) );
              assert(ss_angle >= M_PI/2.0);
                      // We know that this angle must be larger than Pi/2 here by construction of the problem

              length_type scale_center_shell_dist( otherShell_radius * std::sin(M_PI-(angle_diff+ss_angle)) / std::sin(angle_diff) );
              assert(scale_center_shell_dist>0.0);

              // Check whether the calculated distance is within the forseen bounds and warn if not
              if( scale_center_shell_dist > std::sqrt(scale_center_to_shell_y*scale_center_to_shell_y 
                                                      + scale_center_to_shell_z*scale_center_to_shell_z) * (1.0+TOLERANCE) )
              {
        //                 log.warn("Orthogonal cylinder scaling, EDGE_HITS_BARREL case: scale-center-to-shell distance is out of foreseen bounds:");
        //                 log.warn("   distance=%s, scale_center_to_shell_y=%s, scale_center_to_shell_z=%s",
        //                              scale_center_shell_dist, scale_center_to_shell_y, scale_center_to_shell_z);
                      ; // FIXME fix logger first!
              }

              if(scale_angle <= M_PI/4.0)
              {
                  Real cos_scale_angle( std::cos(scale_angle) );
                  z1_new = std::min(z1, (scale_center_z + cos_scale_angle * scale_center_shell_dist));
                  r_new  = std::min(r, (this->*r1_function[di])(z1_new));
              }
              else
              {
                  Real sin_scale_angle( std::sin(scale_angle) );
                  r_new  = std::min(r,  (scale_center_r + sin_scale_angle * scale_center_shell_dist));
                  z1_new = std::min(z1, (this->*z1_function[di])(r_new));
              }
          } 
          break;

          case BARREL_HITS_BARREL:
          {   
              // The scaling cylinder hits the barrel of 'otherShell' with its barrel
              if(scale_angle == 0.0)
              {
                  // In this case the cylinders can never hit; just leave the lenghts as they are
                  r_new = r;
                  z1_new = z1;
              }
              else
              {
                  // case scale_angle > 0
                  r_new  = std::min(r, (ref_to_shell_y - otherShell_radius));
                  z1_new = std::min(z1, (this->*z1_function[di])(r_new));
              }
          }     
          break;

          default:
              
              throw unsupported("get_dr_dzright_dzleft_to_CylindricalShape: Bad situation for making cylinders against cylinders.");
              
      }; // switch(collision_situation)
      
      z2_new = std::min(z2, (this->*z2_function[di])(r_new));
      
      return create_vector<position_type>(r_new, z1_new, z2_new); // TODO what happens with z2_new?
            
    };
    
    /***********************************/
    /* Collision of parallel cylinders */
    /***********************************/
    inline position_type get_dr_dzright_dzleft_to_parallel_CylindricalShape()
    {
        // Calculates the new testShell dimensions when scaled with respect
        // to another shell that is parallel in orientation
        
        length_type r_new(r), z1_new(z1), z2_new(z2); // TODO what happens with z2_new?
        
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
        {   // otherShell collides with the scaling cylinder ('testShell') on top
            z1_new = std::min(z1, (ref_to_shell_z - otherShell_hl) );
            r_new  = std::min(r,  (this->*r1_function[this->di])(z1_new) );
                     // TODO if z1 hasn't changed we also don't have to recalculate this
        }
        else
        {   // otherShell collides with the scaling cylinder ('testShell') on the radial side
            r_new  = std::min(r, (ref_to_shell_r - otherShell_radius) );
            z1_new = std::min(z1, (this->*z1_function[this->di])(r_new) );
        }
        
        return create_vector<position_type>(r_new, z1_new, z2_new); // TODO what happens with z2_new?
    };
    
    /****************************/
    /* FURTHER HELPER FUNCTIONS */
    /****************************/
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
    
    // Functions to feed into the rootfinder
    gsl_function_pt_type edge_hits_edge_h1_eq_pt;
    gsl_function_pt_type edge_hits_edge_r1_eq_pt;
    
    // Important vectors and lengths used throughout the calculations
    // Local coordinate system in orthogonal scaling
    position_type       local_x, local_y, local_z;    
    position_type       ref_to_shell_vec;
    length_type         ref_to_shell_x, ref_to_shell_y, ref_to_shell_z;
    length_type         ref_to_shell_x2, ref_to_shell_y2, ref_to_shell_z2;
    length_type         scale_center_to_shell_x, scale_center_to_shell_y, scale_center_to_shell_z;
    
    length_type         r1_min, h1_min; // used in rootfinder routine in EDGE_HITS_EDGE case
    
    static Logger&      log_;
    
    
  public:       // CONSTRUCTOR
    
    CylinderScalingHelperTools( CylinderScalingFunctionsWrap<Ttraits_> *CSF_, 
                                position_type   testShell_reference_point_,
                                position_type   testShell_orientation_vector_,
                                position_type   testShell_dimensions_,
                                position_type   otherShell_position_t_,
                                position_type   otherShell_orientation_vector_,
                                position_type   otherShell_dimensions_,
                                direction_type  direction_,
                                position_type   scale_center_info_, 
                                position_type   scale_angle_info_                 ):
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
         // Assign the scaling functions to the pointer arrays
         // for the two different directions direction = -1 (di = 0)
         r1_function[0] = &CylinderScalingHelperTools<Ttraits_>::r_left;
         z1_function[0] = &CylinderScalingHelperTools<Ttraits_>::z_left;
         z2_function[0] = &CylinderScalingHelperTools<Ttraits_>::z_right;
         // direction = +1 (di = 1)
         r1_function[1] = &CylinderScalingHelperTools<Ttraits_>::r_right;
         z1_function[1] = &CylinderScalingHelperTools<Ttraits_>::z_right;
         z2_function[1] = &CylinderScalingHelperTools<Ttraits_>::z_left;
         
         // Assign pointers to the functions that are used to solve the implicit
         // equation in the edge-hits-edge intersection case via the GSL rootfinder
         edge_hits_edge_h1_eq_pt = ( double (*)(double, void*) )(&CylinderScalingHelperTools<Ttraits_>::edge_hits_edge_h1_eq);
         edge_hits_edge_r1_eq_pt = ( double (*)(double, void*) )(&CylinderScalingHelperTools<Ttraits_>::edge_hits_edge_h1_eq);
         
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




/*** LOGGER ***/
template<typename Ttraits_>
Logger& CylinderScalingHelperTools<Ttraits_>::log_(Logger::get_logger("ecell.CylinderScalingHelperTools"));






#endif /* CYLINDERSCALING_HPP */
