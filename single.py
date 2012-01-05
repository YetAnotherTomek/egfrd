from _gfrd import (
    SphericalShell,
    CylindricalShell,
    Cylinder,
    Sphere,
    CylindricalSurface,
    PlanarSurface,
    )
from _greens_functions import *
from greens_function_wrapper import *
from constants import EventType
import utils

from domain import (
    Domain,
    ProtectiveDomain)

from shells import *

__all__ = [
    'CylindricalSurfaceSingle',
    'PlanarSurfaceSingle',
    'SphericalSingle',
    'Single',
    'NonInteractionSingle',
    'InteractionSingle',
    'CylindricalSurfaceInteraction',
    'PlanarSurfaceInteraction',
    ]


class Single(ProtectiveDomain):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see 
    coordinate.py. For each coordinate the Green's function is 
    specified.

    """
    def __init__(self, domain_id, shell_id, reactionrules):
        ProtectiveDomain.__init__(self, domain_id, shell_id)

        self.multiplicity = 1                           # all singles have only one particle

        assert self.testShell
        assert isinstance(self.testShell, testSingle)
        self.pid_particle_pair = self.testShell.pid_particle_pair      # the particle in the single
        self.structure = self.testShell.structure                      # the structure in/on which the particle lives

        self.reactionrules = reactionrules              # the reaction rules for mono molecular reactions
        self.rrule = None                               # the reaction rule of the actual reaction (property reactionrule)
        self.k_tot = self.calc_ktot(reactionrules)

    def getD(self):
        return self.pid_particle_pair[1].D
    D = property(getD)

    def get_reaction_rule(self):
        if self.rrule == None:
            self.rrule = self.draw_reaction_rule(self.reactionrules)
        return self.rrule
    reactionrule = property(get_reaction_rule)

    def get_inner_a(self):
        # inner_a is the outer radius of the 'inner shell', which demarkates
        # the space where the CoM of the particle can go. The space can be
        # one, two or three dimensional and can be linear, round, spherical or
        # cylindric in form.
        return self.shell.shape.radius - self.pid_particle_pair[1].radius

    def draw_new_position(self, dt, event_type):
        '''Draw a new position of the particle in the coordinate system given that an
        event of event_type has occured. No periodic boundary conditions are applied yet
        to the coordinates, as the Domains do not understand this.

        Note also that it calculates the positional information PRIOR to any identity or
        structure change.

        When for example a particles dissociates from the membrane (a structure change),
        it is placed next to the membrane. 'draw_new_position' will only determine the
        location in the membrane where it dissociates, not the position where the particle
        should be placed.
        '''
        pass

    def is_reset(self):
        return self.dt == 0.0 and self.event_type == EventType.SINGLE_ESCAPE

    def draw_reaction_time_tuple(self):
        # MOVE to ProtectiveDomain class as this is actually a property of a single particle
        # and not of a Single domain. It should be a method to a particle (or species)
        # but that's too much hassle. Now make it a function (!) of Domain class for
        # easy access. This mechanism is used everywhere a particle can decay.
        """Return a (reaction time, event type)-tuple for a unimolecular reaction.

        """
        if self.k_tot <= 0:
            dt = numpy.inf
        elif self.k_tot == numpy.inf:
            dt = 0.0
        else:
            rnd = myrandom.uniform()
            if rnd == 0:
                dt = numpy.inf
            else:
                dt = (1.0 / self.k_tot) * (- math.log(rnd)) # log(1/x) == - log(x)
        return dt, EventType.SINGLE_REACTION

    def draw_escape_time_tuple(self):
        """Return an (escape time, event type)-tuple.

        In case of an interaction single, handles only the escape in the
        non-interacting coordinate.
        
        """
        if self.getD() == 0:
            dt = numpy.inf
        else:
            dt = draw_time_wrapper(self.greens_function())

        event_type = EventType.SINGLE_ESCAPE    # TODO Event is not an ESCAPE when in 1D there is drift
        return dt, event_type

    def check(self):
    # performs an internal consistency check
        pass

    def __str__(self):
        pid = self.pid_particle_pair[0]
        sid = self.pid_particle_pair[1].sid
        name = self.world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        return 'Single[%s: %s, ST%s, %s]' % (self.domain_id, pid, name,
                                             self.pid_particle_pair[1].position)


class NonInteractionSingle(Single):
    """1 Particle inside a shell, no other particles around. 

    There are 3 types of NonInteractionSingles:
        * SphericalSingle: spherical shell, 3D movement.
        * PlanarSurfaceSingle: cylindrical shell, 2D movement.
        * CylindricalSurfaceSingle: cylindrical shell, 1D movement.

    NonInteractionSingle hold two independent processes:
        * Diffusion of the particle in n (1-3) dimensions
        * Decay process of the particle

    """
    def __init__(self, domain_id, shell_id, reactionrules):
        Single.__init__(self, domain_id, shell_id, reactionrules)

    def initialize(self, t):
        # self.shell.shape.radius = self.pid_particle_pair[1].radius
        self.dt = 0.0
        self.last_time = t
        self.event_type = EventType.SINGLE_ESCAPE


    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple())

    def draw_new_position(self, dt, event_type):
        oldpos = self.pid_particle_pair[1].position

        if self.D == 0:
            newpos = oldpos
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos = oldpos
        else:
            # calculate r
            if event_type == EventType.SINGLE_ESCAPE:
            # Moving this checks to the Green's functions is not a good 
            # idea, because then you'd draw an unused random number.  
            # The same yields for the draw_new_com and draw_new_iv.  

            # TODO Either include check for non-zero drift here, or overload
            # this method in CylindricalSurfaceSingle to check exit interface
                r = self.get_inner_a()
            else:
                gf = self.greens_function()
                r = draw_r_wrapper(gf, dt, self.get_inner_a())
                # Note that in case of 1D diffusion r has a direction. It is
                # the lateral displacement and can be positive or negative. 
                # In other cases r is always positive and denotes
                # radial displacement from the center.

            # calculate other coordinates
            displacement = self.create_position_vector(r)

            # This should be checked in the unit test of random_vector.
            # if __debug__:
            #     scale = self.pid_particle_pair[1].radius
            #     if feq(length(displacement), abs(r), typical=scale) == False:
            #         raise AssertionError('displacement != abs(r): %g != %g.' % 
            #                              (length(displacement), abs(r)))

            # Add displacement to shell.shape.position, not to particle.position.  
            # This distinction is important only in the case of an 
            # asymmetric 1D domain (r0 != 0, or drift), since draw_r always 
            # returns a distance relative to the centre of the shell (r=0), 
            # not relative to r0.

            # note that we need to make sure that the shell.shape.position and displacement vector
            # are in the structure to prevent the particle leaving the structure due to numerical errors
            newpos = self.shell.shape.position + displacement

        return newpos

### TODO move this in some way to the shell making stuff in shells.py
    def calculate_shell_size_to_single(self, closest, distance_to_shell, geometrycontainer):
        # only if the closest shell is also a 'simple' Single
        assert isinstance(closest, NonInteractionSingle)

        min_radius1 = self.get_min_shell_size()

        D1 = self.getD()
        if D1 == 0:
            return min_radius1

        min_radius2 = closest.pid_particle_pair[1].radius
        D2 = closest.getD()
        min_radius12 = min_radius1 + min_radius2
        sqrtD1 = math.sqrt(D1)


        singlepos = self.pid_particle_pair[1].position
        closestpos = closest.shell.shape.position
        distance = geometrycontainer.distance(singlepos, closestpos)

        shell_size = min(sqrtD1 / (sqrtD1 + math.sqrt(D2))
                        * (distance - min_radius12) + min_radius1,
                        distance_to_shell / SAFETY)
        if shell_size < min_radius1:
            shell_size = min_radius1

        return shell_size


### this can probably both go
    def get_max_shell_size(self, geometrycontainer, domains):
        # This calculates the maximum shell size that the single can have in the current
        # situation.
        # TODO This is now general for all NonInteractionSingles and should be different
        # for the different dimensions since they have different domains (spheres vs cylinders)
        # and scale the domains in different directions (r or z).

        # 1. Calculate the maximal dimensions of the shell
        # TODO We now just check in a sphere around the particle and size accordingly. For
        # NonInteractionSingles on the dna and membrane this is quite inefficient-> want to
        # size in the appropriate coordinate (r or z).
        singlepos = self.pid_particle_pair[1].position
#        closest, distance_to_shell = \
#            geometrycontainer.get_closest_obj(singlepos, domains, ignore=[self.domain_id],
#                                              ignores=[self.structure.id])

        neighbor_domains = geometrycontainer.get_neighbor_domains(singlepos, domains, ignore=[self.domain_id, ])
        neighbor_surfaces = geometrycontainer.get_neighbor_surfaces(singlepos, ignores=[self.structure.id])
        neighbors = neighbor_domains + neighbor_surfaces

        # 2. Calculate the maximum allowed radius of the spherical domain
        # check ALL domains, not just nearest neighbor but also next nearest etc, and get the nearest one.
        # In case the neighbor was a NonInteractionSingle we may need to tweak the distance
        distance_domain_closest = (numpy.inf, None)
        for domain, distance in neighbors:
            if isinstance(domain, NonInteractionSingle):
                # leave at least MULTI_SHELL_FACTOR space for another NonInteractionSingle (to make sure that
                # the other single will not have an overlapping multi shell when starting to multi
                distance_max = geometrycontainer.world.distance(singlepos, domain.pid_particle_pair[1].position)
                distance = min(distance, distance_max - domain.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR)
                # TODO check this!
            elif isinstance(domain, PlanarSurface):
                distance = distance + self.pid_particle_pair[1].radius

            # get the minimum of currect and previous minimum
            distance_domain = (distance, domain)
            distance_domain_closest = min(distance_domain_closest, distance_domain)

        distance_to_shell, closest = distance_domain_closest

        # 3. Apply still other rules if the nearest is another NonInteractionSingle (optimization)
        # TODO shell size has no real meaning for cylinderical domains
        if isinstance(closest, NonInteractionSingle):
            new_shell_size = self.calculate_shell_size_to_single(closest, distance_to_shell,
                                                                 geometrycontainer)
        else:  # Pair or Multi or Surface
            new_shell_size = distance_to_shell / SAFETY

        return min(new_shell_size, geometrycontainer.get_max_shell_size())


class SphericalSingle(NonInteractionSingle, NonInteractionSingles, hasSphericalShell):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector:
          theta, phi.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):
        assert isinstance(testShell, SphericalSingletestShell)

        hasSphericalShell.__init__(self, testShell, domain_id)           # here also the shell is created
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def greens_function(self):
        return GreensFunction3DAbsSym(self.getD(),
                                          self.get_inner_a())

    # specific shell methods for the SphericalSingle
    # these return potentially corrected dimensions
    def shell_list_for_single(self):
        min_radius = self.pid_particle_pair[1].radius * MULTI_SHELL_FACTOR
#        min_radius = self.testShell.get_min_radius()
        if  self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            fake_shell = self.create_new_shell(position, min_radius, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_radius = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            fake_shell = self.create_new_shell(position, min_radius, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def create_updated_shell(self, position):
        try:
            radius = self.testShell.determine_possible_shell([self.domain_id], [self.structure.id])
            return self.create_new_shell(position, radius, self.domain_id)

        except ShellmakingError as e:
            raise Exception('SphericalSingle, create_updated_shell failed: %s' %
                            (str(e)))

    def create_position_vector(self, r):
        return random_vector(r)

    def __str__(self):
        return 'Spherical' + Single.__str__(self)


class PlanarSurfaceSingle(NonInteractionSingle, NonInteractionSingles, hasCylindricalShell):
    """1 Particle inside a (cylindrical) shell on a PlanarSurface. (Hockey 
    pucks).

        * Particle coordinates on surface: x, y.
        * Domain: radial r. (determines x and y together with theta).
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        hasCylindricalShell.__init__(self, testShell, domain_id)
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def greens_function(self):
        return GreensFunction2DAbsSym(self.D,
                                          self.get_inner_a())

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    def shell_list_for_single(self):
        # The shell should always be larger that the bare minimum for a test shell
        # Note that in case of a cylindrical shell this fits inside of the spherical multi shell
        min_radius = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
#        min_radius, _, _ = self.testShell.get_min_dr_dzright_dzleft()
        if self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            half_length = self.shell.shape.half_length
            fake_shell = self.create_new_shell(position, min_radius, half_length, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_radius = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.radius < min_radius:
            position = self.shell.shape.position
            half_length = self.shell.shape.half_length
            fake_shell = self.create_new_shell(position, min_radius, half_length, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def create_updated_shell(self, position):
        # TODO what should we do with the position now?
        try:
            dr, dz_right, dz_left = self.testShell.determine_possible_shell([self.domain_id], [self.structure.id])
            center, radius, half_length = self.r_zright_zleft_to_r_center_hl(self.testShell.get_referencepoint(),
                                                                             self.testShell.get_orientation_vector(),
                                                                             dr, dz_right, dz_left)            
            return self.create_new_shell(center, radius, half_length, self.domain_id)
        except ShellmakingError as e:
            raise Exception('PlanarSurfaceSingle, create_updated_shell failed: %s' % str(e) )

    def create_position_vector(self, r):
        # project the vector onto the surface unit vectors to make sure that the coordinates are in the surface
        x, y = random_vector2D(r)
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

    def __str__(self):
        return 'PlanarSurface' + Single.__str__(self)


class CylindricalSurfaceSingle(NonInteractionSingle, NonInteractionSingles, hasCylindricalShell):
    """1 Particle inside a (cylindrical) shell on a CylindricalSurface. 
    (Rods).

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none.
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules):

        hasCylindricalShell.__init__(self, testShell, domain_id)
        NonInteractionSingle.__init__(self, domain_id, shell_id, reactionrules)

    def getv(self):
        return self.pid_particle_pair[1].v
    v = property(getv)

    def get_inner_a(self):
        # Here the free coordinate is not the radius of the domain as in the
        # PlanarSurfaceSingle and SphericalSingle, although it is treated as a
        # radius in NonInteractionSingle.draw_new_position for simplicity
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        inner_half_length = self.get_inner_a()
        return GreensFunction1DAbsAbs(self.D, self.v, 0.0, -inner_half_length, inner_half_length)

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    def shell_list_for_single(self):
        # The shell should always be larger that the bare minimum for a test shell
        # Note that in case of a cylindrical shell this fits inside of the spherical multi shell
        min_half_length = self.pid_particle_pair[1].radius * math.sqrt(MULTI_SHELL_FACTOR**2 - 1.0)
#        _, min_half_length, _ = self.testShell.get_min_dr_dzright_dzleft()
        if self.shell.shape.half_length < min_half_length:
            position = self.shell.shape.position
            radius = self.shell.shape.radius
            fake_shell = self.create_new_shell(position, radius, min_half_length, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_half_length = self.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.half_length < min_half_length:
            position = self.shell.shape.position
            radius = self.shell.shape.radius
            fake_shell = self.create_new_shell(position, radius, min_half_length, self.domain_id)
            return [(self.shell_id, fake_shell), ]
        else:
            return self.shell_list

    def create_updated_shell(self, position):
        # TODO what should we do with the position now?
        # maybe update the reference_point of the shell before updating the shell
        try:
            dr, dz_right, dz_left = self.testShell.determine_possible_shell([self.domain_id], [self.structure.id])
            dz_right = min(dz_right, dz_left)       # make sure the domain is symmetric around the particle
            dz_left  = dz_right                     # This is not necessary but it's assumed later
            center, radius, half_length = self.r_zright_zleft_to_r_center_hl(self.testShell.get_referencepoint(),
                                                                             self.testShell.get_orientation_vector(),
                                                                             dr, dz_right, dz_left)
            return self.create_new_shell(center, radius, half_length, self.domain_id)
        except ShellmakingError as e:
            raise Exception('CylindricalSurfaceSingle, create_updated_shell failed: %s' % str(e) )

    def create_position_vector(self, z):
        if utils.feq(z, self.get_inner_a()):
            # Escape, can be either to the left or to the right.
            # The side of escape should be decided on by the flux through
            # both boundaries at the escape time
            # TODO, include drift
            z = myrandom.choice(-1, 1) * z 

        # project the com onto the surface unit vectors to make sure that the coordinates are in the surface
        return z * self.shell.shape.unit_z

    def __str__(self):
        return 'CylindricalSurface' + Single.__str__(self)


class InteractionSingle(Single, hasCylindricalShell, Others):
    """Interactions singles are used when a particle is close to a 
    surface.

    There are 2 types of InteractionSingles:
        * PlanarSurfaceInteraction.
        * CylindricalSurfaceInteraction.

    The InteractionSingle holds three independent processes:
        * Diffusion in the direction that allows the interaction (radial or z)
        * Diffusion in the direction that is normal the the direction of interaction
        * Decay of the particle into a different or no particle
    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):

        hasCylindricalShell.__init__(self, testShell, domain_id)
        assert self.testShell
        assert isinstance(self.testShell, testInteractionSingle)

        Single.__init__(self, domain_id, shell_id, reactionrules)

        self.interactionrules = interactionrules        # the surface interaction 'reactions'
        self.intrule = None
        self.interaction_ktot = self.calc_ktot(interactionrules)

        self.surface = self.testShell.surface           # The surface with which the particle is
                                                        # trying to interact

    def get_interaction_rule(self):
        if self.intrule == None:
            self.intrule = self.draw_reaction_rule(self.interactionrules)
        return self.intrule
    interactionrule = property(get_interaction_rule)

        # There are defined here but are overloaded later
        # These ProtectiveDomains have a Cylindrical shell that is scaled in r and z
        # Below methods go together with this property
    def get_inner_dz_left(self):
        pass

    def get_inner_dz_right(self):
        pass

    def initialize(self, t):    
        # initialize the domain object with the appropriate time to allow reuse of
        # the Domain at different times
        self.dt = 0
        self.last_time = t
        self.event_type = None

    def determine_next_event(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_time_tuple(),
                   self.draw_reaction_time_tuple(),
                   self.draw_iv_event_time_tuple())     # above two events also occur in NonInteractingSingle
                                                        # IV event is added in InteractionSingle

    def draw_iv_event_time_tuple(self):
        """Return an (interaction/escape time, event type)-tuple.
        Note that an interaction is modelled similarly to an escape.
 
        Note: we are not calling single.draw_iv_event_type() just yet, but 
        postpone it to the very last minute (when this event is executed 
        in fire_single). So IV_EVENT can still be an iv escape or an iv 
        interaction.

        """
        dt = draw_time_wrapper(self.iv_greens_function())

        return dt, EventType.IV_EVENT

    def draw_iv_event_type(self):
        gf = self.iv_greens_function()
        event_kind = draw_event_type_wrapper(gf, self.dt)
        if event_kind == PairEventKind.IV_REACTION:     # Note REACTION -> INTERACTION
            return EventType.IV_INTERACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()

    def __str__(self):
        pass

class PlanarSurfaceInteraction(InteractionSingle):
    """1 Particle close to a PlanarSurface, inside a cylindrical shell 
    placed against the surface with its flat face.

        * Particle coordinates inside shell: r, theta, z.
        * Coordinates: radial r, cartesian z.
        * Initial position: r = 0, z = z0 (relative to center shell).
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):
        # shell_orientation_vector and z0 are supplied by arguments because Domains can't do
        # spatial calculations (they don't know about the world)
        # NOTE shell_orientation_vector should be normalized and pointing from the surface to the particle!

        InteractionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules)

        # z0 is the initial position of the particle relative to the center
        # of the shell in the direction normal to the surface (so z0 is usually negative).
        self.z0 = -self.get_inner_dz_left() + self.testShell.particle_surface_distance

    def get_inner_dz_right(self):
        # This is the distance from the center of the shell
        # to the flat boundary away from the surface
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_dz_left(self):
        # This is the distance from the center of the shell to the surface
        # Note that it is implied that the shell sticks out at the other end
        # of the planar surface by the particle radius, and that the surface
        # has no thickness.
        dist_to_surface = self.shell.shape.half_length - self.pid_particle_pair[1].radius
        return dist_to_surface #- self.pid_particle_pair[1].radius

    def greens_function(self):
        # This is the Greens function that does not produce REACTION event
        #return GreensFunction3DAbsSYm(self.getD(), get_inner_a())
        return GreensFunction2DAbsSym(self.D, self.get_inner_a())

    def iv_greens_function(self):
        # The green's function that also modelles the association of the particle
        # with the planar surface.
        return GreensFunction1DRadAbs(self.D, self.interaction_ktot, self.z0,
                                      -self.get_inner_dz_left(), self.get_inner_dz_right())

    def draw_new_position(self, dt, event_type):
        """Draws a new position for the particle for a given event type at a given time
        Note that it also 'moves the particle into the surface' when an interaction
        takes place.
        """
        oldpos = self.pid_particle_pair[1].position
        if self.D == 0:
            newpos = oldpos
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos = oldpos             # just return oldpos, particle will be removed anyway
        else:
            # calculate r vector
            a = self.get_inner_a()
            if event_type == EventType.SINGLE_ESCAPE:
                r = a
            else:
                gf = self.greens_function()
                r = draw_r_wrapper(gf, dt, a)

            x, y = random_vector2D(r)
            # express the r vector in the unit vectors of the surface to make sure the particle is
            # parallel to the surface (due to numerical problem)
            vector_r = x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

            # calculate z vector
            z_surface = -self.get_inner_dz_left()       # left side of the inner domain
            z_not_surface = self.get_inner_dz_right()   # right side of the inner domain (away from the surface)

            # If the event was not yet fully specified
            if event_type == EventType.IV_EVENT:        
                self.event_type = self.draw_iv_event_type()
                event_type = self.event_type


            if event_type == EventType.IV_INTERACTION:
                z = z_surface 
            elif event_type == EventType.IV_ESCAPE:
                z = z_not_surface
            else:
                gf_iv = self.iv_greens_function()
                z = draw_r_wrapper(gf_iv, dt, z_not_surface, z_surface) # last two args= a, sigma

            # express the vector_z in the unit vectors of the surface to prevent the particle from
            # leaving the surface due to numerical problem
            vector_z = z * self.shell.shape.unit_z

            # The new position is relative to the center of the shell
            # note that we need to make sure that the shell.shape.position and vector_r and vector_z
            # are correct (in the structure) to prevent the particle leaving the structure due to numerical errors
            newpos = self.shell.shape.position + vector_r + vector_z

        return newpos

    def __str__(self):
        return ('PlanarSurfaceInteraction' + Single.__str__(self) + \
               'radius = %s, half_length = %s' %
                (self.shell.shape.radius, self.shell.shape.half_length))


class CylindricalSurfaceInteraction(InteractionSingle):
    """1 Particle close to a CylindricalSurface, inside a cylindrical shell 
    that surrounds the surface.

        * Particle coordinates inside shell: r, theta, z.
        * Domains: composite r-theta, cartesian z.
        * Initial position: r = r, theta = 0, z = z.
        * Selected randomly when drawing displacement vector: none.

    """
    def __init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules):
        InteractionSingle.__init__(self, domain_id, shell_id, testShell, reactionrules, interactionrules)

        # unit_r is the vector from the cylindrical surface to the particle
        # r0 is the distance from the center of the cylinder to the center of the particle
        # z0 is known to be zero (the particle being in the center of the shell in the z direction)
        self.unit_r = self.testShell.reference_vector    
        self.r0     = self.testShell.particle_surface_distance
        self.z0     = 0.0

    def get_inner_dz_left(self):
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_dz_right(self):
        return self.shell.shape.half_length - self.pid_particle_pair[1].radius

    def get_inner_sigma(self):
        return self.surface.shape.radius + self.pid_particle_pair[1].radius

    def greens_function(self):
        # The greens function not used for the interaction but for the other coordinate
        # Free diffusion in z direction, drift is zero. Note that also z0 = 0.0
        return GreensFunction1DAbsAbs(self.D, 0.0, self.z0,
                                      -self.get_inner_dz_left(),
                                      self.get_inner_dz_right())

    def iv_greens_function(self):
        # Green's function used for the interaction
        # Interaction possible in r direction.
        return GreensFunction2DRadAbs(self.D, self.interaction_ktot, self.r0, self.get_inner_sigma(),
                                      self.get_inner_a())

    def draw_new_position(self, dt, event_type):
        oldpos = self.pid_particle_pair[1].position

        if self.D == 0:
            newpos = oldpos
        elif event_type == EventType.SINGLE_REACTION and len(self.reactionrule.products) == 0:
            newpos = oldpos
        else:
            # 1) Draw z.
            if event_type == EventType.SINGLE_ESCAPE:
                # Moving this checks to the Green's functions is not a good 
                # idea, because then you'd draw an unused random number.  
                # The same yields for the draw_new_com and draw_new_iv.  
                # Assume for now that z0 = 0 -> both boundaries are equally likely
                z = myrandom.choice(-1, 1) * self.get_inner_dz_left()
            else:
                gf = self.greens_function()
                z = draw_r_wrapper(gf, dt, self.get_inner_dz_left())

            # Direction matters, so determine the direction of the shell relative to the surface
            # first
#            direction = cmp(numpy.dot(self.shell.shape.unit_z, self.structure.shape.unit_z), 0)
            # express the z_vector into the surface unit vectors to make sure that the coordinates are in the surface
#            z_vector = z * direction * self.structure.shape.unit_z
            z_vector = z * self.shell.shape.unit_z


            # 2) Draw r and theta.
            # If the event was not yet fully specified
            if event_type == EventType.IV_EVENT:        
                self.event_type = self.draw_iv_event_type()
                event_type = self.event_type

            sigma = self.get_inner_sigma()
            a_r = self.get_inner_a()
            gf_iv = self.iv_greens_function()

            if event_type == EventType.IV_ESCAPE:
                r = a_r
            elif event_type == EventType.IV_INTERACTION:
                r = sigma
            else:
                r = draw_r_wrapper(gf_iv, dt, a_r, sigma)
            theta = draw_theta_wrapper(gf_iv, r, dt)

            r_vector = r * rotate_vector(self.unit_r, self.shell.shape.unit_z,
                                     theta)

            # Add displacement to shell.shape.position, not to particle.position.  
            newpos = self.shell.shape.position + z_vector + r_vector

        return newpos

    def __str__(self):
        return ('CylindricalSurfaceInteraction' + Single.__str__(self) + \
               'radius = %s, half_length = %s' %
                (self.shell.shape.radius, self.shell.shape.half_length))

class DummySingle(object):
    def __init__(self):
        self.multiplicity = 1
        self.half_length = 0.0
        self.shellList = [Sphere(NOWHERE, 0.0)]

    def getMinSize(self):
        return 0.0

    def getD(self):
        return 0.0

    def getPos(self):
        return NOWHERE
    pos = property(getPos)

    def __str__(self):
        return 'DummySingle()'

