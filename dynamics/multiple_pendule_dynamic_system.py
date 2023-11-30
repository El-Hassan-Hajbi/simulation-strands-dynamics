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

    def __init__(self, stems, THETA, scheme):
        ## Constructor
        # @param self
        # @param mesh  
        # @param theta0  
        super().__init__()
        self.stems = stems # list of stems : [rods_1st_stem, rods_2nd_stem ...]
        self.THETA = THETA # list of angles for each steam : [THETA_1st_steam, THETA_2nd_steam ...]

        self.scheme = scheme

        # Animations parameters
        self.it = 60.
        self.delta = 0.01
        self.period = 120.
        self.theta0 = THETA[0]

        # Create an array of zeros with the same length
        omega0 = np.zeros_like(self.THETA[0])

        # Concatenate the two arrays horizontally
        self.X = [np.concatenate([omega0, self.THETA[k]]) for k in range(len(self.stems))]

        # Collision detection and response
        self.collisionDetectionOn, self.collisionResponseOn = True, False
    
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
        
        return f
    
    def A_matrix(self, theta):
        N = len(theta)
        A = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if i == j:
                    A[i, j] = 1
                else:
                    A[i, j] = (N-i)/(N-i+1)*np.cos(theta[j] - theta[i])

        return A

    def step(self):
        N = len(self.THETA[0])
        for k in range(len(self.stems)):
            omega, theta = self.X[k][:len(self.X[k])//2], self.X[k][len(self.X[k])//2:]

            A = self.A_matrix(theta)

            f = self.f_vector(theta, omega)
            
            if self.scheme=="explicit":
                self.X[k] += self.delta * F(self.X[k], self.it, A, f)
            elif self.scheme=="implicit":
                self.X[k] = implicit_step(self.X[k], self.delta, self.it, A, f)
            elif self.scheme=="runge-kutta":
                self.X[k] = runge_kutta_step(self.X[k], self.delta, self.it, A, f)
            
            self.THETA[k] = self.X[k][N:]
            
        for j in range(len(self.stems)):
            # Updating positions of all the rods
            self.stems[j][0].positions[2], self.stems[j][0].positions[3] = self.stems[j][0].positions[0] + l*np.sin(self.THETA[j][0]), -l*np.cos(self.THETA[j][0])
            for i, rod in enumerate(self.stems[j][1:]):
                rod.positions[0], rod.positions[1] = self.stems[j][i-1 + 1].positions[2], self.stems[j][i-1 + 1].positions[3] # The first extremite should match with the last extremite of the previous rod
                rod.positions[2], rod.positions[3] = rod.positions[0]+l*np.sin(self.THETA[j][i + 1]), rod.positions[1]-l*np.cos(self.THETA[j][i + 1]) # The other extremite of the ith rod should be updated knowing the updated value of theta
        
        # updating time
        self.it += self.delta

        # Collision
        if self.collisionDetectionOn:
            segs = [[[(self.stems[j][i].positions[0], self.stems[j][i].positions[1]), (self.stems[j][i].positions[2], self.stems[j][i].positions[3])] for i in range(len(self.stems[0]))] for j in range(len(self.stems))]
            collision, intersections_point, stem1IDX, stem2IDX, rod1IDX, rod2IDX = CollisionDetection(segments=segs)
            if collision : 
                print("intersection points = ", intersections_point)
                print("stem 1 IDX of intersection = ", stem1IDX)
                print("stem 2 IDX of intersection = ", stem2IDX)
                print("rod 1 IDX of intersection = ", rod1IDX)
                print("rod 2 IDX of intersection = ", rod2IDX)
                if self.collisionResponseOn:
                # ----------- Compute impulse response
                    Fcol = CollisionResponse()
                    # ----------- Update angular velocities of interesected rods
                    for s in range(len(intersections_point)): 
                        omega1, omega2 = self.X[stem1IDX[s]][:len(self.X[0])//2], self.X[stem2IDX[s]][:len(self.X[0])//2]
                        omega1[rod1IDX[s]] = -omega1[rod1IDX[s]]
                        omega2[rod2IDX[s]] = -omega2[rod2IDX[s]]