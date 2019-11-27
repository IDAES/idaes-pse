
__all__ = ['Parallelepiped','RectPrism','Cube','Rhombohedron','Cuboctahedron']

from functools import partial
from abc import abstractmethod
from copy import deepcopy
import numpy as np
from math import cos,sin,sqrt

from ..util.util import isZero,areEqual
from .transform_func import TransformFunc,ShiftFunc,ScaleFunc,RotateFunc,ReflectFunc

class Shape(object):
    """ """
    DBL_TOL = 1e-5
    DEFAULT_ALIGNMENT = np.array([[1,0,0],
                                  [0,1,0],
                                  [0,0,1]],dtype=float)

    # === STANDARD CONSTRUCTOR
    def __init__(self,Anchor,Alignment=None):
        self._Anchor = Anchor
        self._Alignment = (Shape.DEFAULT_ALIGNMENT if Alignment is None else Alignment)

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        """ """
        if(self.Anchor is None):
            print('alpha')
            return False
        if(type(self.Alignment) is not np.ndarray): 
            print('A')
            return False
        if(self.Alignment.shape != (3,3)): 
            print('B')
            return False
        if(not areEqual(np.linalg.det(self.Alignment),1.0,Shape.DBL_TOL)): 
            print('C')
            return False
        for AlignmentAxis in self.Alignment:
            if(not areEqual(np.linalg.norm(AlignmentAxis),1.0,Shape.DBL_TOL)): 
                print('D')
                return False
        return True

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(type(TransF) is ShiftFunc or
           type(TransF) is ScaleFunc):
            TransF.transform(self._Anchor)
        elif(type(TransF) is RotateFunc or
             type(TransF) is ReflectFunc):
            TransF.transform(self._Anchor)
            TransF.transform(self._Alignment)
            assert(self.isConsistentWithDesign())
        else:
            raise TypeError
        assert(self.isConsistentWithDesign())

    def shift(self,Shift):
        """

        Args:
            Shift: 

        Returns:

        """
        if(type(Shift) is ShiftFunc):
            self.applyTransF(Shift)
        elif(type(Shift) is np.ndarray):
            self.applyTransF(ShiftFunc(Shift))
        else:
            raise TypeError

    def scale(self,Scale,OriginOfScale=None):
        """

        Args:
            Scale: param OriginOfScale:  (Default value = None)
            OriginOfScale:  (Default value = None)

        Returns:

        """
        if(type(Scale) is ScaleFunc):
            self.applyTransF(Scale)
        elif(type(Scale) is np.ndarray):
            self.applyTransF(ScaleFunc(Scale,OriginOfScale))
        elif(type(Scale) is float or
             type(Scale) is int):
            self.applyTransF(ScaleFunc(np.array([Scale,Scale,Scale],dtype=float)))
        else:
            raise TypeError           
       
    def rotate(self,Rotation,OriginOfRotation=None):
        """

        Args:
            Rotation: param OriginOfRotation:  (Default value = None)
            OriginOfRotation:  (Default value = None)

        Returns:

        """
        if(type(Rotation) is RotateFunc):
            self.applyTransF(Rotation)
        elif(type(Rotation) is np.ndarray):
            self.applyTransF(RotateFunc(Rotation,OriginOfRotation))
        else:
            raise TypeError

    def reflect(self,Reflection):
        """

        Args:
            Reflection: 

        Returns:

        """
        if(type(Reflection) is ReflectFunc):
            self.applyTransF(Reflection)
        else:
            raise TypeError
            
    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def isInShape(self,P):
        """

        Args:
            P: 

        Returns:

        """
        raise NotImplementedError

    __contains__ = isInShape

    @abstractmethod
    def getBBox(self):
        """ """
        raise NotImplementedError

    # === BASIC QUERY METHODS
    @property
    def Anchor(self):
        """ """
        return self._Anchor

    @property
    def Alignment(self):
        """ """
        return self._Alignment

class Polyhedron(Shape):
    """ """
    # === STANDARD CONSTRUCTOR
    def __init__(self,V,F,Anchor,Alignment=None):
        Shape.__init__(self,Anchor,Alignment)
        self._V = V
        self._F = F
        self._FacetNorms = self.__calcFacetNorms()
        self._FacetDirections = self.__calcFacetDirections()
        assert(self.isConsistentWithDesign())

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        """ """
        if(type(self.V) is not list): return False
        if(len(self.V) < 4): return False
        if(type(self.F) is not list): return False
        if(len(self.F) < 4): return False
        for Facet in self.F:
            if(type(Facet) is not list): return False
            if(len(Facet) < 3): return False
        return Shape.isConsistentWithDesign(self)

    # === AUXILIARY METHODS
    def __calcFacetNorms(self):
        FacetNorms = []
        for Facet in self.F:
            #import code; code.interact(local=dict(locals(),**globals()));
            Norm = np.cross(self.V[Facet[1]]-self.V[Facet[0]],
                            self.V[Facet[2]]-self.V[Facet[0]])
            Norm /= np.linalg.norm(Norm)
            FacetNorms.append(Norm)
        return FacetNorms

    def __calcFacetDirections(self):
        FacetDirections = []
        for Facet in self.F:
            assert(len(Facet)>=3)
            FacetDirection = np.array([0,0,0],dtype=float)
            for v in Facet:
                FacetDirection += self.V[v]
            FacetDirection /= len(Facet)
            FacetDirection -= self.Anchor
            FacetDirections.append(FacetDirection)
        return FacetDirections

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(type(TransF) is ShiftFunc):
            for v in self._V:
                TransF.transform(v)
        elif(type(TransF) is ScaleFunc): 
            for v in self._V:
                TransF.transform(v)
            for direction in self._FacetDirections:
                TransF.transform(direction)
            for norm in self._FacetNorms:
                norm /= TransF.Scale
        elif(type(TransF) is RotateFunc):
            for v in self._V:
                TransF.transform(v)
            for direction in self._FacetDirections:
                TransF.transformDirection(direction)
            for norm in self._FacetNorms:
                TransF.transformDirection(norm)
        elif(type(TransF) is ReflectFunc):
            raise NotImplementedError #TODO
        else:
            raise TypeError
        Shape.applyTransF(self,TransF)
        assert(self.isConsistentWithDesign())

    # === PROPERTY EVALUATION METHODS
    def isInShape(self,P,tol=Shape.DBL_TOL):
        """

        Args:
            P: param tol:  (Default value = Shape.DBL_TOL)
            tol:  (Default value = Shape.DBL_TOL)

        Returns:

        """
        for f in range(len(self.F)):
            if(not self.satisfiesFacet(P,f,tol)): return False
        return True
    
    __contains__ = isInShape

    def satisfiesFacet(self,P,f,tol=Shape.DBL_TOL):
        """

        Args:
            P: param f:
            tol: Default value = Shape.DBL_TOL)
            f: 

        Returns:

        """
        P01 = P - self.V[self.F[f][0]]
        return (np.inner(P01,self.FacetNorms[f]) < tol)

    def getBounds(self):
        """ """
        MinP = deepcopy(self.V[0])
        MaxP = deepcopy(self.V[0])
        for v in self.V:
            MinP = np.minimum(MinP,v)
            MaxP = np.maximum(MaxP,v)
        return MinP,MaxP

    # === BASIC QUERY METHODS
    @property
    def V(self):
        """ """
        return self._V
    
    @property
    def F(self):
        """ """
        return self._F
    
    @property 
    def FacetNorms(self):
        """ """
        return self._FacetNorms
        
    @property
    def FacetDirections(self):
        """ """
        return self._FacetDirections

class Cuboctahedron(Polyhedron):
    """ """
    # === STANDARD CONSTRUCTOR
    def __init__(self,R,Center=None):
        self._R = sqrt(2)
        V = [np.array([ 1, 0, 1],dtype=float),
             np.array([ 0,-1, 1],dtype=float),
             np.array([-1, 0, 1],dtype=float),
             np.array([ 0, 1, 1],dtype=float),
             np.array([ 1, 1, 0],dtype=float),
             np.array([ 1,-1, 0],dtype=float),
             np.array([-1,-1, 0],dtype=float),
             np.array([-1, 1, 0],dtype=float),
             np.array([ 1, 0,-1],dtype=float),
             np.array([ 0,-1,-1],dtype=float),
             np.array([-1, 0,-1],dtype=float),
             np.array([ 0, 1,-1],dtype=float)]
        F = [[ 0, 3, 2, 1],
             [ 1, 6, 9, 5],
             [ 2, 7,10, 6],
             [ 3, 4,11, 7],
             [ 0, 5, 8, 4],
             [ 8, 9,10,11],
             [ 1, 2, 6],
             [ 2, 3, 7],
             [ 0, 4, 3],
             [ 0, 1, 5],
             [ 6,10, 9],
             [ 7,11,10],
             [ 4, 8,11],
             [ 5, 9, 8]]
        Polyhedron.__init__(self,V,F,np.array([0,0,0],dtype=float))
        self.scale(R/sqrt(2))
        if(Center is not None):
            self.shift(Center)

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign():
        """ """
        #TODO: Actually put something here
        return Polyhedron.isConsistentWithDesign(self)

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(isinstance(TransF,ScaleFunc)):
            if(TransF.isIsometric):
                self._R *= TransF.Scale[0]
            else:
                raise ValueError('Cuboctahedron applyTransF: Can only scale isometrically')
        Polyhedron.applyTransF(self,TransF)
        assert(self.isConsistentWithDesign())


class Parallelepiped(Polyhedron):
    """ """
    # === STANDARD CONSTRUCTOR
    def __init__(self,Vx,Vy,Vz,BotBackLeftCorner=None):
        self._Vx = Vx
        self._Vy = Vy
        self._Vz = Vz
        if(np.inner(np.cross(Vx,Vy),Vz)<Parallelepiped.DBL_TOL):
            raise ValueError("Vx,Vy,Vz must have a positive box product. TODO: Understand why...")
            Vx,Vy = Vy,Vx
        if(BotBackLeftCorner is None):
            BotBackLeftCorner = np.array([0,0,0],dtype=float)
            
        V = [BotBackLeftCorner,
             BotBackLeftCorner+Vx,
             BotBackLeftCorner+Vy,
             BotBackLeftCorner+Vz,
             BotBackLeftCorner+Vx+Vy,
             BotBackLeftCorner+Vx+Vz,
             BotBackLeftCorner+Vy+Vz,
             BotBackLeftCorner+Vx+Vy+Vz]
        F = [[0,1,5,3],
             [0,3,6,2],
             [2,6,7,4],
             [1,4,7,5],
             [3,5,7,6],
             [0,2,4,1]]
        Polyhedron.__init__(self,V,F,Anchor=deepcopy(BotBackLeftCorner))

    # === CONSTRUCTOR - From edge lengths and angles
    @classmethod
    def fromEdgesAndAngles(cls,A,B,C,alpha,beta,gamma,BotBackLeftCorner=None):
        """

        Args:
            A: param B:
            C: param alpha:
            beta: param gamma:
            BotBackLeftCorner: Default value = None)
            B: 
            alpha: 
            gamma: 

        Returns:

        """
        Vx = np.array([A,0,0])
        Vy = np.array([B*np.cos(gamma),B*np.sin(gamma),0])
        zcomp = np.sqrt(1-pow(np.cos(alpha),2)-pow(np.cos(beta),2))
        Vz = np.array([C*np.cos(alpha),C*np.cos(beta),C*zcomp])
        #import code; code.interact(local=dict(locals(),**globals()));
        return cls(Vx,Vy,Vz,BotBackLeftCorner)

    # === CONSTRUCTOR - From POSCAR file header information
    @classmethod 
    def fromPOSCAR(cls,filename):
        """

        Args:
            filename: 

        Returns:

        """
        with open(filename,'r') as infile:
            CommentLine = infile.readline()
            GSLine = infile.readline().split()
            GS = float(GSLine[0])
            VxLine = infile.readline().split()
            Vx = GS*np.array([float(VxLine[0]),
                              float(VxLine[1]),
                              float(VxLine[2])],dtype=float)
            VyLine = infile.readline().split()
            Vy = GS*np.array([float(VyLine[0]),
                              float(VyLine[1]),
                              float(VyLine[2])],dtype=float)
            VzLine = infile.readline().split()
            Vz = GS*np.array([float(VzLine[0]),
                              float(VzLine[1]),
                              float(VzLine[2])],dtype=float)
        return cls(Vx,Vy,Vz)

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithClassDesign(self):
        """ """
        # TODO: Add more checks here
        return Polyhedron.isConsistentWithDesign(self)

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(isinstance(TransF,ScaleFunc)):
            self._Vx *= TransF.Scale[0]
            self._Vy *= TransF.Scale[1]
            self._Vz *= TransF.Scale[2]
        Polyhedron.applyTransF(self,TransF)
        assert(self.isConsistentWithDesign())       

    # === PROPERTY EVALUATION METHODS
    def getFractionalCoords(self,P,blnRelativeToCenter=False):
        """

        Args:
            P: param blnRelativeToCenter:  (Default value = False)
            blnRelativeToCenter:  (Default value = False)

        Returns:

        """
        FracOrigin = (self.getCenter() if blnRelativeToCenter else self.Anchor)
        Omega = self.getVolume()
        fracX = 1.0/Omega*np.inner(P-FracOrigin,np.cross(self.Vy,self.Vz))
        fracY = 1.0/Omega*np.inner(P-FracOrigin,np.cross(self.Vz,self.Vx))
        fracZ = 1.0/Omega*np.inner(P-FracOrigin,np.cross(self.Vx,self.Vy))
        return np.array([fracX,fracY,fracZ])

    def getCenter(self):
        """ """
        return sum(self.V)/len(self.V)

    def getVolume(self):
        """ """
        return np.dot(np.cross(self.Vx,self.Vy),self.Vz)

    # === BASIC QUERY METHODS
    @property
    def Vx(self):
        """ """
        return self._Vx

    @property
    def Vy(self):
        """ """
        return self._Vy

    @property
    def Vz(self):
        """ """
        return self._Vz

    @property
    def isUnit(self):
        """ """
        if(not areEqual(np.linalg.norm(self.Vx),1.0,Shape.DBL_TOL)):
            return False
        if(not areEqual(np.linalg.norm(self.Vy),1.0,Shape.DBL_TOL)):
            return False
        if(not areEqual(np.linalg.norm(self.Vz),1.0,Shape.DBL_TOL)):
            return False
        return True

class Rhombohedron(Parallelepiped):
    """ """
    # === STANDARD CONSTRUCTOR
    def __init__(self,L,Alpha,BotBackLeftCorner=None):
        self._L = L
        self._Alpha = Alpha
        Lx = np.array([L,0,0])
        Ly = np.array([L*cos(Alpha),L*sin(Alpha),0])
        Lz = np.array([L*cos(Alpha),
                       L*(cos(Alpha)-cos(Alpha)**2)/sin(Alpha),
                       L*sqrt(1-3*cos(Alpha)**2+2*cos(Alpha)**3)/sin(Alpha)])
        Parallelepiped.__init__(self,Lx,Ly,Lz,BotBackLeftCorner)

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(isinstance(TransF,ScaleFunc)):
            if(TransF.isIsometric):
                self._L *= TransF.Scale[0]
            else:
                raise ValueError('Rhombohedron applyTransF: Can only scale isometrically')
        Polyhedron.applyTransF(self,TransF)
        assert(self.isConsistentWithDesign())       

    # === PROPERTY EVALUATION METHODS
    @property
    def L(self):
        """ """
        return self._L
        
    @property
    def Alpha(self):
        """ """
        return self._Alpha

class RectPrism(Parallelepiped):
    """ """
    # === STANDARD CONSTRUCTOR
    def __init__(self,Lx,Ly,Lz,BotBackLeftCorner=None):
        self._Lx = Lx
        self._Ly = Ly
        self._Lz = Lz
        Parallelepiped.__init__(self,
                                np.array([Lx, 0, 0],dtype=float),
                                np.array([ 0,Ly, 0],dtype=float),
                                np.array([ 0, 0,Lz],dtype=float),
                                BotBackLeftCorner)

    # === CONSTRUCTOR - Bounding box of a list of points
    @classmethod
    def fromPointsBBox(cls,Pts):
        """

        Args:
            Pts: 

        Returns:

        """
        #import code; code.interact(local=dict(locals(),**globals()));
        minX = min(_[0] for _ in Pts)
        maxX = max(_[0] for _ in Pts)
        minY = min(_[1] for _ in Pts)
        maxY = max(_[1] for _ in Pts)
        minZ = min(_[2] for _ in Pts)
        maxZ = max(_[2] for _ in Pts)
        return cls(maxX-minX,maxY-minY,maxZ-minZ,np.array([minX,minY,minZ]))

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(isinstance(TransF,ScaleFunc)):
            self._Lx *= TransF.Scale[0]
            self._Ly *= TransF.Scale[1]
            self._Lz *= TransF.Scale[2]
        Parallelepiped.applyTransF(self,TransF)
        assert(self.isConsistentWithDesign())       

    # === BASIC QUERY METHODS
    @property
    def Lx(self):
        """ """
        return self._Lx

    @property
    def Ly(self):
        """ """
        return self._Ly

    @property
    def Lz(self):
        """ """
        return self._Lz

class Cube(RectPrism):
    """ """
    # === STANDARD CONSTRUCTOR
    def __init__(self,L,BotBackLeftCorner=None):
        self._L = L
        RectPrism.__init__(self,L,L,L,BotBackLeftCorner)

    # === MANIPULATION METHODS
    def applyTransF(self,TransF):
        """

        Args:
            TransF: 

        Returns:

        """
        if(isinstance(TransF,ScaleFunc)):
            if(TransF.isIsometric):
                self._L *= TransF.Scale[0]
            else:
                raise ValueError('Cube applyTransF: Can only scale isometrically')
        RectPrism.applyTransF(self,TransF)
        assert(self.isConsistentWithDesign())       

    # === BASIC QUERY METHODS
    @property 
    def L(self):
        """ """
        return self._L
