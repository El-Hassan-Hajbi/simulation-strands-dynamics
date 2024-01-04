import numpy as np
from scipy.optimize import fsolve
from .abstract_dynamic_system import AbstractDynamicSystem
from scipy.integrate import odeint 
from collision import CollisionDetection, CollisionResponse

g, l, m, K, rayleigh_coeff, wind_coeff, wind_freq, xref = 9.8, 1, 1, 10, 0.1, 2, 1, 0.1 # N/m

def F(X, t, A, f):
    # X = (w1, .., wN, O1, .., ON)
    H = (np.linalg.inv(A) @ f).flatten()
    N = len(X)//2
    W = X[:N]
    res=np.concatenate([H, W])
    return res

def beta(i, j, N):
    if j < i:
        return N-i+1
    if i == j:
        return 0
    else :
        return N-j+1
        
def alpha(i,j, omega, theta):
    N = len(theta)
    return (l*omega[i]*np.sin(theta[j] - theta[i])*beta(i, j, N))/G(i, theta)
def G(s, theta):
    c = 0
    for j in range(s, len(theta)):
        for i in range(j):
            c+= l*theta[i]*np.sin(theta[i] + theta[s])
        c*= m*l
    return -c    

def implicit_step(X, delta, t, A, f):
    # Define the function to find the new state implicitly
    def equation(Y):
        return X - Y + delta * F(Y, t, A, f)
    
    # Use a numerical solver to find the new state
    X_new = fsolve(equation, X)
    return X_new


def runge_kutta_step(X, delta, t, A, f):
    k1 = F(X, t, A, f)
    k2 = F(X + delta*k1/2, t + delta/2, A, f)
    k3 = F(X + delta*k2/2, t + delta/2, A, f)
    k4 = F(X + delta*k3, t + delta, A, f)

    return X + delta * (k1/6 + k2/3 + k3/3 + k4/6)

## Dummy dynamic system just to test
class MultiplePenduleDynamicSystem(AbstractDynamicSystem):

    def __init__(self, stems, THETA, scheme, viewer):
        ## Constructor
        # @param self
        # @param mesh  
        # @param theta0  
        super().__init__()
        self.stems = stems # list of stems : [rods_1st_stem, rods_2nd_stem ...]
        self.THETA = THETA # list of angles for each steam : [THETA_1st_steam, THETA_2nd_steam ...]

        self.scheme = scheme
        self.viewer = viewer

        # Animations parameters
        self.it = 60.
        self.delta = 0.005
        self.period = 120.
        self.theta0 = THETA[0]
        self.cpt = 0
        self.e = 1 # coefficient of elasticity / restitution
        # Create an array of zeros with the same length
        omega0 = np.zeros_like(self.THETA[0])

        # Concatenate the two arrays horizontally
        self.X = [np.concatenate([omega0, self.THETA[k]]) for k in range(len(self.stems))]

        # Collision detection and response
        self.collisionDetectionOn, self.collisionResponseOn = True, True
    
    def f_vector(self, theta, omega):
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
            wind_action = wind_coeff*np.sin(wind_freq*self.it)*s*(l*np.sin(theta[s]) - xref)
            c = (c + tmp + ressort_elastic + rayleigh_damping + wind_action) # to add ressort  + 2*k*theta[s]/(m*l**2)
            f[s] = -c
        
        return m*l**2*f
    
    def A_matrix(self, theta):
        N = len(theta)
        A = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if i == j:
                    A[i, j] = m*l**2
                else:
                    A[i, j] = m*l**2*(N-i)/(N-i+1)*np.cos(theta[j] - theta[i])

        return A

    def step(self):
        # Degrees of freedom
        N = len(self.THETA[0])
        # Collision
        fc = [N*[0] for k in range(len(self.stems))]
        #print(len(fc[0]))
        if self.collisionDetectionOn:
            
            segs = [[[(self.stems[j][i].positions[0], self.stems[j][i].positions[1]), (self.stems[j][i].positions[2], self.stems[j][i].positions[3])] for i in range(len(self.stems[0]))] for j in range(len(self.stems))]
            boolean, collisions, index_maping = CollisionDetection(segments=segs, viewer=self.viewer)
            if boolean :     
                self.cpt+=1
                print(f"The {self.cpt}th collision")
                for col in collisions:
                    col.print()
                if self.collisionResponseOn:
                # ----------- Compute impulse response
                    rc = CollisionResponse(ndl=N, collisions=collisions, index_mapping=index_maping, THETA=self.THETA, dt=self.delta, qqDot=self.X)
                    f_collision = self.e*rc/self.delta # rc = dt * fc
                    for i, index_solid_contact in enumerate(list(index_maping)):
                        fc[index_solid_contact] = f_collision[i*N:N*(i+1)]
                        
        #print(len(fc[0]))

        # Dynamics + rc (Signorini - Coulomb)          
        for k in range(len(self.stems)):
            omega, theta = self.X[k][:len(self.X[k])//2], self.X[k][len(self.X[k])//2:]

            M = self.A_matrix(theta)

            f = self.f_vector(theta, omega)
            
            #print(k, len(f), len(fc[k]))
            # scheme
            if self.scheme=="explicit":
                self.X[k] += self.delta * F(self.X[k], self.it, M, f+fc[k])
            elif self.scheme=="implicit":
                self.X[k] = implicit_step(self.X[k], self.delta, self.it, M, f+fc[k])
            elif self.scheme=="runge-kutta":
                self.X[k] = runge_kutta_step(self.X[k], self.delta, self.it, M, f+fc[k])
            
            self.THETA[k] = self.X[k][N:]
            
        # Update positions
        for j in range(len(self.stems)):
            # Updating positions of all the rods
            self.stems[j][0].positions[2], self.stems[j][0].positions[3] = self.stems[j][0].positions[0] + l*np.sin(self.THETA[j][0]), -l*np.cos(self.THETA[j][0])
            for i, rod in enumerate(self.stems[j][1:]):
                rod.positions[0], rod.positions[1] = self.stems[j][i-1 + 1].positions[2], self.stems[j][i-1 + 1].positions[3] # The first extremite should match with the last extremite of the previous rod
                rod.positions[2], rod.positions[3] = rod.positions[0]+l*np.sin(self.THETA[j][i + 1]), rod.positions[1]-l*np.cos(self.THETA[j][i + 1]) # The other extremite of the ith rod should be updated knowing the updated value of theta
        
        # updating time
        self.it += self.delta

          