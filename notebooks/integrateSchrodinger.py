
import numpy as np
from scipy.fftpack import fft,ifft
from scipy.integrate import solve_ivp

class IntegrateSchrodinger:
    '''
    Class to integrate the Schrodinger equation for arbitrary 1D potentials.
    We use the Kosloff method based on the FFT
    '''
    def __init__(self, x, psi0, V, k0 = None, hbar=1, m=1, t0=0.0, dt=0.01):
        
        # Validation of array inputs
        self.x, psi0, self.V = map(np.asarray, (x, psi0, V))
        self.N = self.x.shape[0]
        assert self.x.shape == (self.N,)
        assert psi0.shape == (self.N,)
        assert self.V.shape == (self.N,)
        
        # Set internal parameters
        self.hbar = hbar
        self.m = m
        self.t = t0
        self.dt = dt
        self.dx = self.x[1] - self.x[0]
        self.dk = 2 * np.pi / (self.N * self.dx)
        
        # Set minimum momentum
        if k0 == None:
            self.k0 = -np.pi/self.dx
        else:
            self.k0 = k0
        # Set momentum grid
        self.k = self.k0 + self.dk * np.arange(self.N)
        
        # Define psi_x and psi_x_old
        self.psi_x = psi0
        self.psi_x_old = psi0 # Store old value of psi_x for finite difference
        self.psi_x_all = psi0.reshape(1,-1)
        self.dt = dt
        

    
    def compute_psi_k_from_psi_x(self, psi_x):
        '''
        FFT to psi(x,t)
        We define psi(x,t) multiplied by the term so that psiF_x and psiF_k are Fourier pairs
        '''
        psiF_x = (psi_x * np.exp(-1j * self.k[0] * self.x) * self.dx / np.sqrt(2 * np.pi))
        
        psiF_k = fft(psiF_x)
        
        psi_k = psiF_k * np.exp(-1j * self.x[0] * self.dk * np.arange(self.N))
        
        return psi_k


    def compute_psi_x_from_psi_k(self, psi_k):
        '''
        Inverse FFT to psi(k,t)
        We define psi(x,t) multiplied by the term so that psiF_x and psiF_k are Fourier pairs
        '''
        psiF_k = psi_k* np.exp(1j * self.x[0] * self.dk * np.arange(self.N))
        
        psiF_x = ifft(psiF_k)
        
        psi_x = (psiF_x* np.exp(1j * self.k[0] * self.x) * np.sqrt(2 * np.pi) / self.dx)
        
        return psi_x
    
    def diff_RK(self, t, psi_x):
        '''
        Function to pass to the Runge-Kutta method
        '''
        # 1. Calculate nabla^2psi(x,t) by:
        # a) Obtain psi(k,t) using FFT
        psi_k = self.compute_psi_k_from_psi_x(psi_x)

        # b) Multiply psi(k,t) by hbar^2 k^2/(2m)
        nabla_psi_k = -psi_k*self.k**2

        # c) Obtain the kinetic part in the x space using the inverse FFT
        nabla_psi_x = -self.hbar**2/(2*self.m)*self.compute_psi_x_from_psi_k(nabla_psi_k)

        # 2. Add V(x)psi(x,t) to obtain H*psi(x,t)
        H_psi = nabla_psi_x +  self.V*self.psi_x
        
        return -1j/self.hbar*H_psi

        
    def evolve(self, t_increase, store_all = False):
        '''
        Algorithm to integrate the Schrodinger equation up to time t_final
        '''
        if t_increase<0:
            print('Error. Time increment must be positive')
            return
        
        N_steps = int(np.floor((t_increase)/self.dt)+1) #Number of steps to reach t_final from current time
        
        #Initialization: Obtain psi_1 from psi_0 using a Runge-Kutta method
        sol = solve_ivp(self.diff_RK, [0, self.dt], self.psi_x_old)
        self.psi_x = sol.y[:,-1]
        if store_all:
            self.psi_x_all = np.concatenate((self.psi_x_all, self.psi_x.reshape(1,-1)), axis=0)
    
        #For each time step:
        for i in range(1, N_steps):
            print('Time step: {}/{}'.format(i, N_steps), end='\r')
        # 1. Calculate nabla^2psi(x,t) by:
            # a) Obtain psi(k,t) using FFT
            psi_k = self.compute_psi_k_from_psi_x(self.psi_x.copy())
            
            # b) Multiply psi(k,t) by hbar^2 k^2/(2m)
            nabla_psi_k = -psi_k*self.k**2
            
            # c) Obtain the kinetic part in the x space using the inverse FFT
            nabla_psi_x = -self.hbar**2/(2*self.m)*self.compute_psi_x_from_psi_k(nabla_psi_k)

        # 2. Add V(x)psi(x,t) to obtain H*psi(x,t)
            H_psi = nabla_psi_x +  self.V*self.psi_x
            
                
        # 3. Calculate psi(x, t+dt) = psi(x, t-dt) - 2i/h*dt*H*psi(x,t) 
            aux = self.psi_x.copy() # Store value to set it as the old value later
            self.psi_x = self.psi_x_old - 2*1j/self.hbar*self.dt*H_psi
            self.psi_x_old = aux #Set old value of psi
           
            if store_all:
                self.psi_x_all = np.concatenate((self.psi_x_all, self.psi_x.reshape(1,-1)), axis=0)
        
        # Increase final time
        self.t += t_increase        



def calc_norm(phi, xmin=-5, xmax=5, n_points=100):
    '''
    Calculates the norm of an eigenfunction
    Args:
        phi (np.array): wave function values
        xmin (float): minimum value of x
        xmax (float): maximum value of x
        n_points (int): Number of grid points of x
    Returns:
        (float): norm of phi(x)
    '''
    h = (xmax - xmin)/n_points
    return 1./np.sqrt(np.sum(np.conjugate(phi)*phi*h))
     

def empirical_energy1D( phi, potential, xmin=-8, xmax = 8,
                       n_points=200, hbar=1, m=1):
    '''
    Function to calculate the empirical energy of a wavefunction
    Args:
      phi (np.array): Wavefunctions
      potential (np.array): potential V(x)
      xmin (int): minimum value of x
      xmax (int): maximum value of x
      n_points (int): number of grid points
      hbar (float): h bar
      m (float): mass
    Returns:
      E (np.array): empirical energies
    '''
    # Normalize phi just in case
    h = (xmax - xmin)/n_points

    def energy(phi,potential,h):
        '''
        Calculates the empirical energy for one wavefunction
        Args:
        phi (np.array): Wavefunctions
        potential (np.array): potential V(x)
        h (float): lattice size
        Returns:
        E (float): empirical energy 
        '''
        C = 1./np.sqrt(np.sum(np.conjugate(phi)*phi*h))
        phi = C*phi
        # We first calculate the second derivative of phi
        phir = np.concatenate(( phi[1:], np.zeros(1)), axis=0) # We add 0 at the extrema. It makes sense because phi(x)->0 at x->+-inf
        phil = np.concatenate(( np.zeros(1), phi[:-1]), axis=0)

        deriv = (phir - 2*phi + phil)/(h*h)
        return np.sum((-hbar*hbar/(2*m)*np.conjugate(phi)*deriv + potential*(np.conjugate(phi)*phi))*h)

    E = np.array([energy(phi[i,:], potential[i,:], h) for i in range(phi.shape[0])])   
    return E