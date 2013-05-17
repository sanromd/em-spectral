#!/usr/bin/env python
# encoding: utf-8
#
import numpy as np

class Dimension(object):
    r"""
    basic class to represent a dimension
    """
    
    # ========== Property Definitions ========================================
    @property
    def delta(self):
        r"""(float) - distance between grid points"""
        return (self.upper-self.lower) / float(self.num_points)
    @property
    def grid_points(self):
        r"""(ndarrary(:)) - Location of grid points"""
        if self._grid_points is None:
            self._grid_points = np.linspace(self.lower,self.upper,self.num_points+1)
        return self._grid_points
    _grid_points = None
    @property
    def centers(self):
        r"""(ndarrary(:)) - Location of all cell center coordinates
        for this dimension"""
        if self._centers is None:
            self._centers = np.linspace(self.lower+self.delta/2,self.upper-self.delta/2,self.num_points)
        return self._centers
    _centers = None

    def __init__(self,params,dim=0):
        r"""
        Creates a Dimension object
        
        See :class:`Dimension` for full documentation
        """

        self.name = params.dimension_name[dim]
        self.num_points = params.dimension_points[dim]
        self.lower = params.dimension_lower[dim]
        self.upper= params.dimension_upper[dim]
        self.num_cells = self.num_points - 1

    def __str__(self):
        output = "Dimension %s" % self.name
        output += ":  (num_cells,delta,[lower,upper]) = (%s,%s,[%s,%s])" \
            % (self.num_points,self.delta,self.lower,self.upper)
        return output

class Simple_dimension(object):
    r"""
    basic class to represent a dimension
    
    x = geometry.dimension(name,x_lower,x_upper,n)

    where:

     - *name*   - (string) string Name of dimension
     - *lower*  - (float) Lower extent of dimension
     - *upper*  - (float) Upper extent of dimension
     - *n*      - (int) Number of points

    Example:

    >>> from geometry import dimension
    >>> x = dimension('x',0.,1.,100)
    >>> print x
    Dimension x:  (num_cells,delta,[lower,upper]) = (100,0.01,[0.0,1.0])
    >>> x.name
    'x'
    >>> x.num_points
    100
    >>> x.delta
    0.01
    >>> x.grid_points[0]
    0.0
    >>> x.grid_points[1]
    0.01
    >>> x.grid_points[-1]
    1.0
    >>> x.centers[-1]
    0.995
    >>> len(x.centers)
    100
    >>> len(x.edges)
    101
    """
    
    # ========== Property Definitions ========================================
    @property
    def delta(self):
        r"""(float) - distance between grid points"""
        return (self.upper-self.lower) / float(self.num_points)
    @property
    def grid_points(self):
        r"""(ndarrary(:)) - Location of grid points"""
        if self._grid_points is None:
            self._grid_points = np.linspace(self.lower,self.upper,self.num_points+1)
        return self._grid_points
    _grid_points = None
    @property
    def centers(self):
        r"""(ndarrary(:)) - Location of all cell center coordinates
        for this dimension"""
        if self._centers is None:
            self._centers = np.linspace(self.lower+self.delta/2,self.upper-self.delta/2,self.num_points)
        return self._centers
    _centers = None

    def __init__(self, *args, **kargs):
        r"""
        Creates a Dimension object
        
        See :class:`Dimension` for full documentation
        """
        
        # ========== Class Data Attributes ===================================
        self.name = 'x'
        r"""(string) Name of this coordinate dimension (e.g. 'x')"""
        self.num_points = None
        r"""(int) - Number of cells in this dimension :attr:`units`"""
        self.lower = 0.0
        r"""(float) - Lower computational dimension extent"""
        self.upper = 1.0
        r"""(float) - Upper computational dimension extent"""
        
        # Parse args
        if isinstance(args[0],float):
            self.lower = float(args[0])
            self.upper = float(args[1])
            self.num_points = int(args[2]+1)
            self.num_cells = int(args[2])
        elif isinstance(args[0],basestring):
            self.name = args[0]
            self.lower = float(args[1])
            self.upper = float(args[2])
            self.num_points = int(args[3]+1)
            self.num_cells = int(args[3])
        else:
            raise Exception("Invalid initializer for dimension.")
        
        for (k,v) in kargs.iteritems():
            setattr(self,k,v)

    def __str__(self):
        output = "Dimension %s" % self.name
        if self.units:
            output += " (%s)" % self.units
        output += ":  (num_cells,delta,[lower,upper]) = (%s,%s,[%s,%s])" \
            % (self.num_points,self.delta,self.lower,self.upper)
        return output

class Grid(object):
    r"""
    Basic representation of a single grid in
    
    :Dimension information:
    
        Each dimension has an associated name with it that can be accessed via
        that name such as ``grid.x.num_cells`` which would access the x dimension's
        number of cells.
    
    :Properties:

        If the requested property has multiple values, a list will be returned
        with the corresponding property belonging to the dimensions in order.
         
    :Initialization:
    
        Input:
         - *dimensions* - (list of :class:`Dimension`) Dimensions that are to 
           be associated with this grid
            
        Output:
         - (:class:`grid`) Initialized grid object
    """

    # ========== Property Definitions ========================================
    @property
    def num_dim(self):
        r"""(int) - Number of dimensions"""
        return len(self._dimensions)
    @property
    def dimensions(self):
        r"""(list) - List of :class:`Dimension` objects defining the 
                grid's extent and resolution"""
        return [getattr(self,name) for name in self._dimensions]
    @property
    def num_cells(self): 
        r"""(list) - List of the number of cells in each dimension"""
        return self.get_dim_attribute('num_cells')
    @property
    def lower(self):
        r"""(list) - Lower coordinate extents of each dimension"""
        return self.get_dim_attribute('lower')
    @property
    def upper(self):
        r"""(list) - Upper coordinate extends of each dimension"""
        return self.get_dim_attribute('upper')
    @property
    def delta(self):
        r"""(list) - List of computational cell widths"""
        return self.get_dim_attribute('delta')
    @property
    def grid_points(self):
        r"""(list) - List of grid_points"""
        return self.get_dim_attribute('grid_points')
    @property
    def centers(self):
        r"""(list) - List of center coordinate arrays"""
        return self.get_dim_attribute('centers')
    @property
    def num_points(self):
        "number of points in each dimension"
        return self.get_dim_attribute('num_points')

       
    
    # ========== Class Methods ===============================================
    def __init__(self,dimensions,method='curvilinear',Lnpi=4):
        # Dimension parsing
        dime = dimensions
        if isinstance(dimensions,Dimension):
            dimensions = [dimensions]
        self._dimensions = []
        for dim in dimensions:
            self.add_dimension(dim)

        if method=='spectral':
            print self._dimensions
            self.cart_to_spectral(dime,npis=Lnpi)
        #super(Grid,self).__init__()
    
    # ========== Dimension Manipulation ======================================
    def add_dimension(self,dimension):
        r"""
        Add the specified dimension to this patch
        
        :Input:
         - *dimension* - (:class:`Dimension`) Dimension to be added
        """

        # Add dimension to name list and as an attribute
        if dimension.name in self._dimensions:
            raise Exception('Unable to add dimension. A dimension'\
             +' of the same name: {name}, already exists.'\
             .format(name=dimension.name))

        self._dimensions.append(dimension.name)
        setattr(self,dimension.name,dimension)
        
        
    def get_dim_attribute(self,attr):
        r"""
        Returns a tuple of all dimensions' attribute attr
        """
        return [getattr(getattr(self,name),attr) for name in self._dimensions]
    
    
    # ========== Copy functionality ==========================================
    def __copy__(self):
        return self.__class__(self)
        
    # ========== Grid Operations =============================================
    def curvilinear(self):
        r"""
        grid transformation in curvilinear coordinates
        """
        pass
    def cart_to_curv(self,dimensions):
        r"""
        Apply transformation from cartesian to curvilinar coordinates as specificied in function curvilinear
        """
        pass
    def curv_to_cart(self,dimensions):
        r"""
        back transformation from curvilinear to cartesian
        """
        pass
    def cart_to_spectral(self,dimensions,npis=4):
        r"""
        transformation to a spectral grid [-2pi, 2pi]
        """
        from geometry import Simple_dimension as _dime
        print len(self._dimensions)
        if len(self._dimensions)>=2:
            self.cart_length = self._spectral_off = self._spectral_ratio = self._spectral_n_box = np.zeros([len(self._dimensions)])
            self._spectral_lower = self._spectral_upper = np.zeros([len(self._dimensions)])
            for i in range(0,len(self._dimensions)):
                self.cart_length[i] = (dimensions[i].upper - dimensions[i].lower) 
                self._spectral_ratio[i] = npis*np.pi/self.cart_length[i]
                self._spectral_off[i] = -self.cart_length[i]/2.0 + np.abs(dimensions[i].lower)
                self._spectral_n_box[i] = np.floor(self.cart_length[i]/(2.0*np.pi))
                if dimensions[i].upper>0:
                    self._spectral_lower[i] = (dimensions[i].lower - self._spectral_off[i])*self._spectral_ratio[i]
                    self._spectral_upper[i] = (dimensions[i].upper - self._spectral_off[i])*self._spectral_ratio[i]
                else:
                    self._spectral_lower[i] = (dimensions[i].lower + self._spectral_off[i])*self._spectral_ratio[i]
                    self._spectral_upper[i] = (dimensions[i].upper + self._spectral_off[i])*self._spectral_ratio[i]

                ei_temp = _dime(dimensions[i].name + '_spectral', self._spectral_lower[i],self._spectral_upper[i],dimensions[i].num_points)
                self.add_dimension(ei_temp)
        else:
            print dimensions.name
            self.cart_length = (dimensions.upper - dimensions.lower) 
            self._spectral_ratio = npis*np.pi/self.cart_length
            self._spectral_off = -self.cart_length/2.0 + np.abs(dimensions.lower)
            self._spectral_n_box = np.floor(self.cart_length/(2.0*np.pi))
            if dimensions.upper>0:
                self._spectral_lower = (dimensions.lower + self._spectral_off)*self._spectral_ratio
                self._spectral_upper = (dimensions.upper + self._spectral_off)*self._spectral_ratio
            else:
                self._spectral_lower = (dimensions.lower - self._spectral_off)*self._spectral_ratio
                self._spectral_upper = (dimensions.upper - self._spectral_off)*self._spectral_ratio

            ei_temp = _dime(dimensions.name + '_spectral', self._spectral_lower,self._spectral_upper,dimensions.num_points)
            self.add_dimension(ei_temp)
