import IMP.core
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.algebra
import IMP.em
import IMP.isd


class VesicleMembraneRestraint(IMP.pmi.restraints.RestraintBase):
    """
    renamed from the original restraint name as MembraneRestraint
    IMP.pmi.restraints.basic.MembraneRestraint(representation, objects_inside=inside, objects_above=above,
                                                thickness=40)
    """
    def __init__(self,
                 hier,
                 objects_above=None,
                 objects_inside=None,
                 objects_below=None,
                 center = 0.0,
                 thickness=30.0,
                 softness=3.0,
                 plateau=0.0000000001,
                 resolution=1,
                 weight = 1.0,
                 radius = 1000,
                 label = None):

        """ Setup Membrane restraint

        Simple sigmoid score calculated for particles above,

        @param objects_inside list or tuples of objects in membrane (e.g. ['p1', (10, 30,'p2')])
        @param objects_above list or tuples of objects above membrane
        @param objects_below list or tuples of objects below membrane
        @param thickness Thickness of the membrane along the z-axis
        @param radius Radius of the vesicle
        @param softness Softness of the limiter in the sigmoid function
        @param plateau Parameter to set the probability (=1- plateau)) at the plateau
                       phase of the sigmoid
        @param weight Weight of restraint
        @param label A unique label to be used in outputs and
                     particle/restraint names.
        input a list of particles, the slope and theta of the sigmoid potential
        theta is the cutoff distance for a protein-protein contact
        """

        self.hier = hier
        model = self.hier.get_model()

        rname = "MembraneRestraint"
        super(VesicleMembraneRestraint, self).__init__(
            model, name="MembraneRestraint", label=label, weight=weight)

        self.center = center
        self.thickness = thickness
        self.radius = radius
        self.softness = softness
        self.plateau = plateau
        self.linear = 0.02
        self.resolution = resolution

        # Create nuisance particle
        p = IMP.Particle(model)
        z_center = IMP.isd.Nuisance.setup_particle(p)
        z_center.set_nuisance(self.center)

        # Setup restraint
        mr = IMP.pmi.MembraneRestraint(model,
                                       z_center.get_particle_index(),
                                       self.thickness,
                                       self.softness,
                                       self.plateau,
                                       self.linear)

        # Particles above
        if objects_above:
            for obj in objects_above:
                if isinstance(obj, tuple):
                    self.particles_above = self._select_from_tuple(obj)

                elif isinstance(obj, str):
                    self.particles_above = self._select_from_string(obj)
                mr.add_particles_above(self.particles_above)

        # Particles inside
        if objects_inside:
            for obj in objects_inside:
                print "What are objects_inside: ", obj
                if isinstance(obj, tuple):
                    self.particles_inside = self._select_from_tuple(obj)
                    print "Particles inside", self.particles_inside

                elif isinstance(obj, str):
                    self.particles_inside = self._select_from_string(obj)
                    print "vesicle inside", self.particles_inside

                mr.add_particles_inside(self.particles_inside)


        # Particles below
        if objects_below:
            for obj in objects_below:
                if isinstance(obj, tuple):
                    self.particles_below = self._select_from_tuple(obj)

                elif isinstance(obj, str):
                    self.particles_below = self._select_from_string(obj)
                mr.add_particles_below(self.particles_below)

        self.rs.add_restraint(mr)

    def get_particles_above(self):
        return self.particles_above

    def get_particles_inside(self):
        return self.particles_inside

    def get_particles_below(self):
        return self.particles_below

    def _select_from_tuple(self, obj):
        particles = IMP.atom.Selection(self.hier,
                                       molecule = obj[2],
                                       residue_indexes = range(obj[0], obj[1]+1, 1),
                                       resolution = self.resolution).get_selected_particles()

        return particles

    def _select_from_string(self, obj):
        particles = IMP.atom.Selection(self.hier,
                                       molecule = obj,
                                       resolution = self.resolution).get_selected_particles()
        return particles

    # Related to Visualization purpose only
    # It is not straightforward to generate visualization files

    def create_membrane_density(self, file_out='membrane_localization.mrc'):

        """
        Just for visualization of spherical vesicles.
        Writes density of the membrane as mrc files.
        """

        offset = 20.0 * self.thickness
        apix = 10.0
        resolution = 50.0

        # Create a density header of the requested size
        bbox = IMP.algebra.BoundingBox3D(
            IMP.algebra.Vector3D(-self.center - offset, -self.center - offset, -self.center - offset, ),
            IMP.algebra.Vector3D(self.center + offset, self.center + offset, self.center + offset))

        dheader = IMP.em.create_density_header(bbox, apix)
        dheader.set_resolution(resolution)
        dmap = IMP.em.SampledDensityMap(dheader)

        for vox in range(dmap.get_header().get_number_of_voxels()):
            c = dmap.get_location_by_voxel(vox)
            if self._is_membrane(c[2]) == 1:
                dmap.set_value(c[0], c[1], c[2], 1.0)
            else:
                dmap.set_value(c[0], c[1], c[2], 0.0)

        IMP.em.write_map(dmap, file_out)

    def create_vesicle_membrane_density(self, file_out='vesicle_membrane_out.mrc'):

        '''
        Just for visualization purposes.
        Writes density of membrane localization
        '''

        offset = 5.0 * self.thickness
        apix = 3.0
        resolution = 5.0

        # Create a density header of the requested size
        bbox = IMP.algebra.Sphere3D(
            IMP.algebra.Vector3D(self.center + offset, self.center + offset, self.center + offset),
            self.radius)

        dheader = IMP.em.create_density_header(bbox, apix)
        dheader.set_resolution(resolution)
        dmap = IMP.em.SampledDensityMap(dheader)

        for vox in range(dmap.get_header().get_number_of_voxels()):
            c = dmap.get_location_by_voxel(vox)
            if self._is_membrane(c[2]) == 1:
                dmap.set_value(c[0], c[1], c[2], 1.0)
            else:
                dmap.set_value(c[0], c[1], c[2], 0.0)

        IMP.em.write_map(dmap, file_out)

    def _is_membrane(self, z):
        if (z-self.center) < self.thickness/2.0 and  (z-self.center) >= -self.thickness/2.0 :
            return 1
        else:
            return 0


class COMDistanceRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self,
                 root_hier,
                 protein0,
                 protein1,
                 distance=60.0,
                 strength=1.0,
                 label=None,
                 weight=1.0):

        """ Setup an upper-bound distance restraint
            between the center-of-mass of two proteins
        @ param root_hier
        @ param protein1
        @ param protein2
        @ param distance
        @ param strength
        @ param label
        @ param weight
        """

        self.root_hier = root_hier

        model = self.root_hier.get_model()

        rname = "COMDistanceRestraint"
        super(COMDistanceRestraint, self).__init__(
            model, name="COMDistanceRestraint", label=label, weight=weight)

        self.strength = strength
        self.weight = weight
        self.distance = distance

        # Setup restraint
        self.rs = self._create_restraint_set()

        s0 = IMP.atom.Selection(root_hier, molecule= protein0)

        p0 = s0.get_selected_particles()
        if len(p0) == 0:
            print("COMDistanceRestraint: WARNING> cannot select protein %s)" % (prot0))
            exit()
        s1 = IMP.atom.Selection(root_hier, molecule=protein1)
        p1 = s1.get_selected_particles()
        if len(p1) == 0:
            print("COMDistanceRestraint: WARNING> cannot select protein %s)" % (prot0))
            exit()
        # Get COMs
        self.com0 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(model), p0)
        self.com1 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(model), p1)

        coor0 = IMP.core.XYZ(self.com0).get_coordinates()
        coor1 = IMP.core.XYZ(self.com1).get_coordinates()

        d = IMP.algebra.get_distance(coor0, coor1)

        # Distance restraint
        hub = IMP.core.HarmonicUpperBound(self.distance, self.strength)

        df = IMP.core.DistancePairScore(hub)
        dr = IMP.core.PairRestraint(model, df, (self.com0, self.com1))
        self.rs.add_restraint(dr)

    def get_output(self):

        output = super(COMDistanceRestraint, self).get_output()
        return output
