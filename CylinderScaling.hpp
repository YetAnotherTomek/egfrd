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
  
  enum collision_type { BARREL_HITS_FLAT, EDGE_HITS_EDGE, BARREL_HITS_EDGE, 
                        FLAT_HITS_BARREL, EDGE_HITS_BARREL, BARREL_HITS_BARREL };
  
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
    
    // Helper method to determine how precisely the testShell hits the otherShell
    inline collision_type determine_collision_type()
    {
      
      int quadrant = 0;
      collision_type collision_situation;
      
      if(ref_to_shell_x2 < 0.0 && ref_to_shell_y2 < 0.0)
      {
            // Quadrant 1
            quadrant = 1;

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
            quadrant = 2;

            length_type scale_center_to_shell_x_minus_half_length( scale_center_to_shell_x - otherShell_hl );
                        // just because this is reused several times below

            // The case scale_angle=0, i.e. radius remaining constant at scaling, has to be treated separately
            // because in this case the mathematics in the standard case misdetect the collision situation
            if(scale_angle == 0.0)
            {                
                if( std::sqrt( scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length
                               + scale_center_to_shell_y*scale_center_to_shell_y) < r )
                // The lowest point of the static cylinder is within the radius/flat side circle of the scaling cylinder.
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
                    r1_min = (scale_center_to_shell_x - otherShell_hl)*(1.0+TOLERANCE);
                    h1_min = r1_min/tan_scale_angle;
                }
            }
      }
      else if(ref_to_shell_x2 < 0.0 && ref_to_shell_y2 >= 0.0)
      {
            // Quadrant 3
            quadrant = 3;

            // The case scale_angle=0, i.e. radius remaining constant at scaling, has to be treated separately
            // because in this case the mathematics in the standard case misdetect the collision situation
            if(scale_angle == 0.0)
            {
                if(scale_center_to_shell_y < r)
                // The lowest point of the static cylinder is within the radius/flat side circle of the scaling cylinder.
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
                                                            // This is dangerous, permanently changing the property!
                
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
            quadrant = 4;

            assert(ref_to_shell_x2 >= 0.0 && ref_to_shell_y2 >= 0.0);

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
                    r1_min = std::sqrt(scale_center_to_shell_x_minus_half_length*scale_center_to_shell_x_minus_half_length +
                                       scale_center_to_shell_y_minus_radius*scale_center_to_shell_y_minus_radius) * (1.0+TOLERANCE);
                    h1_min = r1_min / tan_scale_angle;
                }
            }
      }
      
      return collision_situation;
      
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
    
    length_type         r1_min, h1_min; // used in rootfinder routine in EDGE_HITS_EDGE case
    
    
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
