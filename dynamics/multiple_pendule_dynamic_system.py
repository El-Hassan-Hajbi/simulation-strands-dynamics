import numpy as np
from scipy.optimize import fsolve
from .abstract_dynamic_system import AbstractDynamicSystem
from scipy.integrate import odeint 

g, l, m, k = 9.8, 0.2, 1, 10 # N/m

def F(X, t, A, f):
    # X = (w1, .., wN, O1, .., ON)
    H = (np.linalg.inv(A) @ f).flatten()
    N = len(X)//2
    #print("N = ", N)
    W = X[:N]
    #print("H = ",H)
    #print("W = ",W)
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
    return l*omega[i]*np.sin(theta[j] - theta[i])*beta(i, j, N)
def implicit_step(X, delta, t, A, f):
    # Define the function to find the new state implicitly
    def equation(Y):
        return X - Y + delta * F(Y, t, A, f)
    
    # Use a numerical solver to find the new state
    X_new = fsolve(equation, X)
    return X_new

## Dummy dynamic system just to test
class MultiplePenduleDynamicSystem(AbstractDynamicSystem):

    def __init__(self, rods, THETA, scheme):
        ## Constructor
        # @param self
        # @param mesh  
        # @param theta0  
        super().__init__()
        self.rods = rods
        self.THETA = THETA
        self.THETA_deriv = [0]*len(THETA) # vecteur derive des angles
        self.THETA_acc = [0]*len(THETA) # vecteur acceleration des angles

        self.scheme = scheme
        # Animations parameters
        self.it = 60.
        self.delta = 0.01
        self.period = 120.
        self.theta0 = THETA[0]
       # Create an array of zeros with the same length
        omega0 = np.zeros_like(self.THETA)

        # Concatenate the two arrays horizontally
        self.X = np.concatenate([omega0, self.THETA])

    
    def step(self):
        
        N = len(self.THETA)
        A = np.zeros((N, N))
        omega, theta = self.X[:len(self.X)//2], self.X[len(self.X)//2:]
        for i in range(N):
            for j in range(N):
                if i == j:
                    A[i, j] = 1
                else:
                    A[i, j] = alpha(i=i,j=j, omega=omega, theta=theta)

        f = np.zeros_like(self.THETA)
        for s in range(N):
            c = 0
            for j in range(s, N):
                c+=g*l*omega[s]*np.sin(theta[s])
                for i in range(j):
                    c+=l*l*omega[i]*np.cos(theta[i]-theta[s])-2*l*omega[s]*omega[i]*np.sin(theta[i])*np.sin(theta[s])
            f[s] = m*c - 2*k*theta[s]
        #A = np.random.rand(N, N)  # Remplacez cela par votre matrice A
        #f = np.random.rand(N, 1)  # Remplacez cela par votre vecteur f
        if self.scheme=="explicit":
            self.X += self.delta * F(self.X, self.it)
        elif self.scheme=="implicit":
            self.X = implicit_step(self.X, self.delta, self.it, A, f)
        
        self.THETA = self.X[len(self.THETA):]
        
        #print("dynamic =", self.THETA)
        
        # Updating positions of all the rods
        self.rods[0].positions[2], self.rods[0].positions[3] = l*np.sin(self.THETA[0]), -l*np.cos(self.THETA[0])
        for i, rod in enumerate(self.rods[1:]):
            rod.positions[0], rod.positions[1] = self.rods[i-1 + 1].positions[2], self.rods[i-1 + 1].positions[3] # The first extremite should match with the last extremite of the previous rod
            rod.positions[2], rod.positions[3] = rod.positions[0]+l*np.sin(self.THETA[i + 1]), rod.positions[1]-l*np.cos(self.THETA[i + 1]) # The other extremite of the ith rod should be updated knowing the updated value of theta
        
        # updating time
        self.it += self.delta



        """
        # Updating params theta using lagrangian + system resolution
        def Lagrangian(theta, theta_deriv, theta_acc):
        
            Params 
            @ theta : vecteur des angles
            @ theta_deriv : vecteur derive des angles
            @ theta_acc : vecteur acceleration des angles
       
            v = [0]*len(theta) # vecteur norme au carré des vitesse generalisee 
            a = v # vecteur derive de la norme au carré des vitesse generalisee 

            v[0] = (l*theta_deriv[0])**2
            for i in range(1, len(theta)):
                v[i] = v[i-1] + (l*theta_deriv[i])**2 + 2*np.cos(theta[i] - theta[i-1])*np.sqrt(v[i-1])*l*theta_deriv[i]
                a[i] = a[i-1] + 2*l**2*theta_deriv[i]*theta_acc[i] - 2*(theta_deriv[i] - theta_deriv[i-1])*np.sin(theta[i] - theta[i-1])*l*theta_deriv[i]*v[i-1] + 2*np.cos(theta[i] - theta[i-1])*l*v[i-1]*theta_acc[i] + 2*np.cos(theta[i] - theta[i-1])*l*a[i-1]*theta_deriv[i]

    def system_pendulum(self, u, t, m, l, g, k):
        #du = derivatives
        # u = variables
        # p = parameters
        # t = time variable
        
        du = np.zeros(2*len(self.THETA))
    
        
        c = np.cos(u[:-2:2] - u[2::2]) # intermediate variables cos(thetai - theta(i+1))
        s = np.sin(u[:-2:2] - u[2::2])  # intermediate variables sin(thetai - theta(i+1))

        
        du[0:-2:2] = u[1:-1:2]   # d(theta) = omega
        du[1] = (2*k/(m*l**2))*u[0] - (g/l)*np.sin(u[0])# dd(theta0)
        for i in range(1, len(self.THETA)):
            du[2*i+1] = ( m2*g*np.sin(u[2])*c - m2*s*(L1*c*u[1]**2 + L2*u[3]**2) - (m1+m2)*g*np.sin(u[0]) ) /( L1 *(m1+m2*s**2) ) # acceleration angulaire
        
        return du       
        """