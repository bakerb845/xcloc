#!/usr/bin/env python3
"""
Purpose: Computes Green's functions in an acoustic homogeneous medium
         assuming an infinite line source in z.
Copyright: Ben Baker distributed under the MIT license.
"""
from math import pow
from math import pi
from numpy import sqrt
from numpy import copy
from numpy import linspace
from numpy import exp
from numpy import asarray
from numpy.linalg import norm
from numpy import ones
from numpy import zeros
from numpy.fft import rfft
from numpy.fft import irfft
from numpy.fft import rfftfreq
import matplotlib.pyplot as plt

class Acoustic2D:
    def __init__(self, vel=5000.0, rho=2700., Q=9999.):
        self.setModelProperties(vel, rho, Q)
        return


    def setModelProperties(self, vel=5000.0, rho=2700., Q=9999.):
        """!
        @brief Sets the model properties.
        @param vel   Model velocity (m/s).
        @param rho   Model density (kg/m**3).
        @param Q     Quality factor.
        @retval 0 indicates success.
        """
        if (vel <= 0.0 or rho <= 0.0 or Q <= 0):
            if (vel <= 0.0):
                print("vel=%f must be positive"%vel)
            if (rho <= 0.0):
                print("rho=%f must be positive"%rho)
            if (Q < 1.0):
                print("Q should be greater than 1"%Q)
            return -1
        self.vel = vel
        self.rho = rho
        self.Q = Q 
        return 0

    def computeGreensLineSource(self, npts, dt, xs, xr, stf=None):
        """!
        @brief Computes the Green's functions in a homogeneous acoustic
               medium assuming.
               .
        @param npts   Number of points in Green's functions.
        @param dt     Sampling period (s) of Green's functions.
        @param xs     The (x,y,z) of the source.
        @param xr     A [nrec x 3] matrix of (x,y,z) positions of receivers.
        @param stf    The source time function.  If not None then this must
                      have dimension [npts].
        @retval G     On success this is a matrix of dimension [nrec x npts]
                      containing the Green's function for each receiver.
        @retval G     On failure this is None.
        """
        G = None
        nrec = 1
        if (xr.shape[0] > 1):
            nrec = xr.shape[0]
        else:
            self.setSourceAndReceiver(xs, xr)
        omega = rfftfreq(npts, d=dt)*2.0*pi # Omega in angular frequency
        nomega = len(omega)
        if (stf is None):
            stfF = ones(nomega, dtype='complex128') # FT of a dirac delta
        else:
            if (npts != len(stf)):
                print("Error npts=%d != len(stf) = %d"%(npts, len(stf)))
                return G
            npts = len(stf)
            stfF = rfft(stf)
        # Set 0 and Nyquist to be 0
        stfF[0] = complex(0.0, 0.0)
        stfF[nomega-1] = complex(0.0, 0.0)
        G = zeros([nrec, npts])
        for irec in range(nrec):
            if (nrec > 1):
                self.setSourceAndReceiver(xs, xr[irec,:])
            Gfd = self.computeLineSourceFD(omega, stfF)
            G[irec,:] = irfft(Gfd, npts)
        return G 

    def setSourceAndReceiver(self, xs, xr):
        """!
        @brief Sets the source and receiver position for this calculation.
        @param xs  The (x,y,z) position of the source.
        @param xr  The (x,y,z) posiiton of the receiver.
        """
        if (len(xs) != 3 or len(xr) != 3):
            if (len(xs) != 3):
                print("Length of xs must be 3")
            if (len(xr) != 3):
                print("Length of xr must be 3")
            return -1
        self.xs = copy(xs)
        self.xr = copy(xr)
        return 0
 
    def computeLineSourceFD(self, omega, stf=None):
        """!
        @brief Computes the analytic response for a infinite line source in
               z in a 2D medium in x and y in the frequency domain.
        @param omega  Angular frequencies at which to evaluate Green's 
                      functions.
        @param stf    The source time function evaluated omega.
        @retval G     The frequency domain Green's functions evaluated at
                      omega and convolved with the source time function stf.
        """
        vel = self.vel
        rho = self.rho
        Q = self.Q
        xr = self.xr
        xs = self.xs
        expipi4 = complex(sqrt(2.0)/2.0, sqrt(2.0)/2.0)
        # Compute the distance - the distance in z is excluded b/c the line
        # source extends infinitely in z-
        r = sqrt( pow(xs[0] - xr[0], 2) + pow(xs[1] - xr[1], 2) )
        coeff = complex(0.0, -1.0/(4.0*rho*vel*vel))*expipi4
        G = zeros(len(omega), dtype='complex')
        for i in range(len(omega)):
            xscal = 0.0
            if (abs(omega[i]) > 0):
                xscal = sqrt(2.0*vel/(pi*omega[i]*r))
            omegar = omega[i]*r
            carg  = complex(-omegar/(2.0*vel*Q), -omegar/vel)
            G[i] = coeff*(xscal*exp(carg))
        # Loop
        if (not stf is None):
            G = stf*G 
        return G

    def ricker(self, npts, dt, peakFreq, lnorm=True, lshift=True):
        """!
        @brief Computes a Ricker wavelet.
        @param npts    Number of points in trace.
        @param dt      Sampling period (s)
        @param peakFreq  Dominant frequency of wavelet (Hz).
        @param lshift    If true then the wavelet is shifted to the start
                         of the trace.
        @param lnorm     If true then the wavelet is normalized by the
                         energy (wavelet/sqrt(dot(wavelet, wavelet))
        @retval ricker   On succesful exit this is a Ricker wavelet stored
                         stored in a array of dimension [npts].
        """
        ricker = zeros(npts)
        pi2f2 = pow(pi*peakFreq, 2)
        npts2 = npts/2
        t = linspace(-npts2*dt, npts2*dt, npts)
        t2 = t*t
        pi2f2t2 = pi2f2*t2
        ricker = (1.0 - 2.0*pi2f2t2)*exp(-pi2f2t2)
        xmax = max(abs(ricker))
        if (lshift):
            tol = 1.0 - 0.995
            work = zeros(npts) 
            for i in range(len(ricker)):
                if (abs(ricker[i]) > tol*xmax): 
                    ncopy = npts - i
                    work[0:ncopy] = ricker[i:npts]
                    break
            ricker = copy(work)
        if (lnorm):
            area = norm(ricker)
            ricker = ricker/area 
        return ricker 
 
if __name__ == "__main__":
    dt = 1.0/6000.0 # Sampling period
    tmodel = 0.5  # Modeling time in seconds /dt + 1)
    vel = 3100.0  # Velocity
    rho = 2700.0  # Density
    Q = 900.0     # Qp
    xscal = 1.e10 # Scalar moment
    npts = int(tmodel/dt + 1) # Length of signals
    ac2d = Acoustic2D(vel, rho, Q)
    
    xs = asarray([0, 0, 0])
    xr = asarray([[100, 100, 0], [200, 100, 0]])
    t = linspace(0, tmodel, npts)
    ricker = ac2d.ricker(npts, dt, 50.0)
    ricker = ricker*xscal
    #plt.plot(t, ricker)
    #plt.show()
    G = ac2d.computeGreensLineSource(npts, dt, xs, xr, stf=ricker)
    #plt.plot(t, G[0,:])
    #plt.plot(t, G[1,:])
    #plt.show()
