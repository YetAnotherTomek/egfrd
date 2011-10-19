#ifndef CUBOIDAL_REGION_HPP
#define CUBOIDAL_REGION_HPP

#include <boost/bind.hpp>
#include "Region.hpp"
#include "Box.hpp"
#include "freeFunctions.hpp"

template<typename Ttraits_>
class CuboidalRegion
    : public BasicRegionImpl<Ttraits_, Box<typename Ttraits_::world_type::traits_type::length_type> >
{
public:
    typedef BasicRegionImpl<Ttraits_, Box<typename Ttraits_::world_type::traits_type::length_type> > base_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::shape_type shape_type;
    typedef typename base_type::rng_type rng_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;
    typedef typename Ttraits_::world_type::species_type species_type;
    typedef std::pair<position_type, position_type> position_pair_type;

    identifier_type const& id() const
    {
        return base_type::id_;
    }

    virtual position_type random_position(rng_type& rng) const
    {
        return ::random_position(base_type::shape(),
                boost::bind(&rng_type::uniform, rng, -1., 1.));
    }

    virtual position_type random_vector(length_type const& r, rng_type& rng) const
    {
        return normalize(
            create_vector<position_type>(
                rng.uniform(-1., 1.),
                rng.uniform(-1., 1.),
                rng.uniform(-1., 1.)), r);
    }

    virtual position_type bd_displacement(length_type const& mean, length_type const& r, rng_type& rng) const
    {
        return create_vector<position_type>(
            rng.normal(0., r),
            rng.normal(0., r),
            rng.normal(0., r));
    }

    virtual length_type drawR_gbd(Real const& rnd, length_type const& r01, Real const& dt, 
                                    Real const& D01, Real const& v) const
    {
         return drawR_gbd_3D(rnd, r01, dt, D01);
    }

    virtual Real p_acceptance(Real const& k_a, Real const& dt, length_type const& r01, position_type const& ipv, 
                                Real const& D0, Real const& D1, Real const& v0, Real const& v1) const
    {
         return k_a * dt / ((I_bd_3D(r01, dt, D0) + I_bd_3D(r01, dt, D1)) * 4.0 * M_PI);
    }

    virtual position_type dissociation_vector( rng_type& rng, length_type const& r01, Real const& dt, 
                                                Real const& D01, Real const& v ) const
    {
        return random_vector( drawR_gbd(rng.uniform(0., 1.), r01, dt, D01, v ), rng );
    }

    virtual Real reaction_volume( length_type const& r0, length_type const& r1, length_type const& rl ) const
    {
        length_type const r01( r0 + r1 );
        length_type const r01l( r01 + rl );
        length_type const r01l_cb( r01l * r01l * r01l );
        length_type const r01_cb( r01 * r01 * r01 );

        return 4.0/3.0 * M_PI * ( r01l_cb - r01_cb );
    }

    virtual position_type surface_dissociation_vector( rng_type& rng, length_type const& r0, length_type const& rl ) const
    {
        return position_type(); //No surface dissociation 'from' the bulk
    }

    virtual position_pair_type geminate_dissociation_positions( rng_type& rng, species_type const& s0, species_type const& s1, position_type const& op, 
        length_type const& rl ) const
    {
        length_type const r01( s0.radius() + s1.radius() );
        Real const D01( s0.D() + s1.D() );

        Real const X( rng.uniform(0.,1.) );

        length_type const r01l( r01 + rl );
        length_type const r01l_cb( r01l * r01l * r01l );
        length_type const r01_cb( r01 * r01 * r01 );
        
        length_type const diss_vec_length( cbrt( X * (r01l_cb - r01_cb ) + r01_cb ) );

        position_type const m( random_vector( diss_vec_length, rng ) );
        
        return position_pair_type( op - m * s0.D() / D01,
                                    op + m * s1.D() / D01 );
    }
    
    virtual position_pair_type special_geminate_dissociation_positions( rng_type& rng, species_type const& s_surf, species_type const& s_bulk, 
        position_type const& op_surf, length_type const& rl ) const
    {
        return position_pair_type(); //No special geminate dissociation 'from' the bulk
    }

    virtual void accept(ImmutativeStructureVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeStructureVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    CuboidalRegion(identifier_type const& id, shape_type const& shape)
        : base_type(id, shape) {}
};

#endif /* CUBOIDAL_REGION_HPP */
