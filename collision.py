import numpy as np
from pydfcp.src.fischerburmeister import *
from pydfcp.src.nsnewton import *
from utils import H_c, M_c
import scipy

def calculer_direction(x1, y1, x2, y2):
    direction = np.array([x2 - x1, y2 - y1])
    return direction / np.linalg.norm(direction)

def calculer_normale(direction):
    # Tourner la direction de 90 degrés (rotation dans le sens horaire)
    normale = np.array([direction[1], -direction[0]])
    return normale / np.linalg.norm(normale)

def find_intersection(segment1, segment2):
    x1, y1 = segment1[0]
    x2, y2 = segment1[1]
    x3, y3 = segment2[0]
    x4, y4 = segment2[1]

    # Parametric equations for the lines
    denum = ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1))
    if denum != 0:
        ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denum
        ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denum

        # Check if the intersection point is within the parameter range of both segments
        if 0 < ua < 1 and 0 < ub < 1:
            # Calculate the intersection point
            intersection_x = x1 + ua * (x2 - x1)
            intersection_y = y1 + ua * (y2 - y1)
            return True, (intersection_x, intersection_y)
    
    # Segments do not intersect
    return False, None

class Collision:
    def __init__(self, t, n, solids, rods):
        self.t = t
        self.n = n
        self.solids = solids
        self.rods = rods

def CollisionDetection(segments):
    """
    @params:
    - segments: a list of segments, each segment is defined as two points (its upper and lower bound)
    
    Returns:
    - (collision: boolean, collisions: list of Collision objects, index_mapping: dictionary)
    """
    collisions = []
    index_mapping = {}
    
    for u in range(len(segments) - 1):
        for v in range(u + 1, len(segments)):
            for i in range(len(segments[u])):
                for j in range(len(segments[v])):
                    collision, intersection_point = find_intersection(segments[u][i], segments[v][j])
                    if collision:
                        t = np.array([1.0, 0.0])  # Replace with your actual tangent vector
                        n = np.array([0.0, 1.0])  # Replace with your actual normal vector
                        #seg = segments[u][i]
                        #t = calculer_direction(*seg[0], *seg[1]) # Replace with your actual tangent vector
                        #n = calculer_normale(t)
                        collision_obj = Collision(t, n, (u, v), (i, j))
                        collisions.append(collision_obj)

                        # Assign normalized indices
                        if u not in index_mapping:
                            index_mapping[u] = len(index_mapping)
                        if v not in index_mapping:
                            index_mapping[v] = len(index_mapping)

    return len(collisions) != 0, collisions, index_mapping
g, l, m, K, rayleigh_coeff, wind_coeff, wind_freq, xref = 9.8, 1, 1, 10, 0.1, 2, 1, 0.1 # N/m

def f_vector(theta, omega, dt):
        N = len(theta)
        f = np.zeros_like(theta)
        for s in range(N):
            c = (g/l)*np.sin(theta[s])
            tmp = 0
            rayleigh_damping = 0
            for j in range(N):
                if j != s:
                    tmp += omega[j]**2*np.sin(theta[j] - theta[s])
                # -------- Rayleigh Damping ---------
                if j <= s:
                    rayleigh_damping -= -rayleigh_coeff*l*N*np.cos(theta[j] - theta[s])*omega[j]
                if j > s:
                    rayleigh_damping -= -rayleigh_coeff*l*(N - j + s)*np.cos(theta[j] - theta[s])*omega[j]
            tmp *= (N-s)/(N-s+1)
            ressort_elastic = 2*K*(theta[s] - np.pi)/(m*l**2)
            wind_action = wind_coeff*np.sin(wind_freq*dt)*s*(l*np.sin(theta[s]) - xref)
            c = (c + tmp + ressort_elastic + rayleigh_damping + wind_action) # to add ressort  + 2*k*theta[s]/(m*l**2)
            f[s] = -c
        
        return m*l**2*f
from scipy.sparse import csc_matrix, csr_matrix
def CollisionResponse(ndl, collisions, index_mapping, THETA,dt,qqDot, verbose=0):
    if(verbose):
        print("Number of collisions: ",len(collisions))
        print("Number of solids in collisions: ",len(index_mapping))
        print("nDoF: ",len(THETA[0]))
    
    Hc = H_c(nc=len(collisions), ndl=ndl, nsc=len(index_mapping), collisions=collisions, index_mapping=index_mapping, THETA=THETA)
    Mc = M_c(ndl=ndl, nsc=len(index_mapping), index_mapping=index_mapping, THETA=THETA)
    if(verbose):
        print("Shape of Hc: ",Hc.shape)
        print("Shape of Mc: ",Mc.shape, type(Mc))
    A = W = Hc @ scipy.sparse.linalg.inv(csc_matrix(Mc)) @ Hc.T
    # On parcours tous les solides en colision 
    Fext, v = np.zeros(ndl*len(index_mapping)), np.zeros(ndl*len(index_mapping))
    for k, solid_index in enumerate(list(index_mapping)):
        omega, theta = qqDot[k][:len(qqDot[k])//2], qqDot[k][len(qqDot[k])//2:] # qqDot du solid k
        Fext[k:ndl+k] = f_vector(omega,theta,dt)
        v[k:ndl+k] = omega

    b = Hc @ (v + scipy.sparse.linalg.inv(csc_matrix(Mc))@Fext*dt)  

    if(verbose):
        print("b :", b)

    # Instantiate FischerBurmeister class
    nContacts = len(collisions)  # Adjust as needed
    mu = 0.5  # Adjust as needed
    A = A.toarray()
    if(verbose):
        print("type of A: ", type(A))
        print("type of b: ", type(b))
    fischer_burmeister = FischerBurmeister(nContacts, mu, A, b)

    # Define the objective function for optimization
    F = fischer_burmeister

    # Initial guess for the solution
    x0 = np.zeros(A.shape[1])

    # Instantiate NonSmoothNewton solver
    solver = NonSmoothNewton(tolerance=1e-6, maxIter=100, verbose=False)

    # Solve the optimization problem
    rc, _ = solver.solve(F, x0)

    Rc = Hc.T @ rc
    if(verbose):
        print("dual shape of rc: ", Rc.shape)
        print("dual of rc: ", Rc)

    return Rc

if __name__ == "__main__":
    # Example usage:
    segment1 = [(0, 0), (2, 2)]
    segment2 = [(1, 1), (3, 3)]
    segment3 = [(0, 2), (2, 0)]

    collision, intersection_point = CollisionDetection([segment1, segment3])

    if collision:
        print("Segments intersect at:", intersection_point)
    else:
        print("Segments do not intersect.")