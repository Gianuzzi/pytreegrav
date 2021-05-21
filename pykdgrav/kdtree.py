<<<<<<< HEAD
from numba import int32, deferred_type, optional, float64, boolean, int64, njit, jit, prange, types, int32
=======
from numba import int32, deferred_type, optional, float64, boolean, int64, njit, jit, prange, types
>>>>>>> e8e00f1d143f6996d6a385dbeb00019c6daa8d01
from numba.experimental import jitclass
import numpy as np
from numpy import empty, empty_like, zeros, zeros_like, sqrt
from numba.typed import List
import gc

node_type = deferred_type()

spec = [
    ('bounds', float64[:,:]), # shape (3,2), contains the coordinate boundaries of the box of the node
    ('size', float64), # largest dimension of the node
    ('delta', float64), # distance from the geometric center of the node to the center of mass of material in the node
    ('indices',int32[:]), # indices of the points located within the node - stored temporarily during treebuild, then deallocated
#    ('points', float64[:,:]),
#    ('masses', float64[:]),
    ('Npoints', int64), # 
    ('hmax', float64),
#    ('softening', float64[:]),
    ('mass', float64),
    ('COM', float64[:]),
    ('IsLeaf', boolean),
    ('HasLeft', boolean),
    ('HasRight', boolean),
    ('left', optional(node_type)),
    ('right', optional(node_type)),
]

@jitclass(spec)
class KDNode(object):
    def __init__(self, idx, x, m, h):
        self.bounds = empty((3,2)) 
        self.indices = idx
        points = x[idx] #np.take(x, idx, axis=0)
        masses = m[idx] #np.take(m, idx)
        softening = h[idx] #np.take(h, idx)

        self.bounds[0,0] = points[:,0].min()
        self.bounds[0,1] = points[:,0].max()
        self.bounds[1,0] = points[:,1].min()
        self.bounds[1,1] = points[:,1].max()
        self.bounds[2,0] = points[:,2].min()
        self.bounds[2,1] = points[:,2].max()
        
#        self.softening = softening
        self.hmax = softening.max()
        
        self.size = max(self.bounds[0,1]-self.bounds[0,0],self.bounds[1,1]-self.bounds[1,0],self.bounds[2,1]-self.bounds[2,0])
#        self.points = points
        self.Npoints = points.shape[0]
#        self.masses = masses
        self.mass = np.sum(masses)
        self.delta = 0.
        if self.Npoints == 1:
            self.IsLeaf = True
            self.COM = points[0]
        else:
            self.IsLeaf = False
            self.COM = zeros(3)
            for k in range(3):
                for i in range(self.Npoints):
                    self.COM[k] += points[i,k]*masses[i]
                self.COM[k] /= self.mass
                self.delta += (0.5*(self.bounds[k,1]+self.bounds[k,0]) - self.COM[k])**2
            self.delta = sqrt(self.delta)

        self.HasLeft = False
        self.HasRight = False        
        self.left = None
        self.right = None

    def GenerateChildren(self, axis, x, m, h):
        if self.IsLeaf:
            return 0
        
        x_axis = x[self.indices,axis]
        med = (self.bounds[axis,0] + self.bounds[axis,1])/2
        index = (x_axis<med)
        if np.any(index):
            self.left = KDNode(self.indices[index], x, m, h)
            self.HasLeft = True
        index = np.invert(index)
        if np.any(index):
            self.right = KDNode(self.indices[index], x, m, h)
            self.HasRight = True
        self.indices = empty(0,dtype=np.int32)
        return 1

node_type.define(KDNode.class_type.instance_type)

@njit
def ConstructKDTree(x, m, softening):
    if len(np.unique(x[:,0])) < len(x):
        raise Exception("Non-unique particle positions are currently not supported by the tree-building algorithm. Consider perturbing your positions with a bit of noise if you really want to proceed.")
    idx = np.arange(len(x), dtype=np.int32)
    root = KDNode(idx, x, m, softening)
    nodes = [root,]
    axis = 0
    divisible_nodes = 1
    count = 0
    while divisible_nodes > 0:
        N = len(nodes)
        divisible_nodes = 0
        for i in range(count, N): # loop through the nodes we spawned in the previous pass
            count += 1
            if nodes[i].IsLeaf:
                continue                
            else:
                generated_children = nodes[i].GenerateChildren(axis, x, m, softening)
                divisible_nodes += generated_children
                if nodes[i].HasLeft:
                    nodes.append(nodes[i].left)
                if nodes[i].HasRight:
                    nodes.append(nodes[i].right)
                    
        axis = (axis+1)%3
    return root
            
