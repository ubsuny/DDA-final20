import numpy as np
import scipy.constants
def getAjk(rk,rj,k):
    """Returns element A_jk of the A matrix
        Args:
           rj = 1x3 numpy array of position of dipole
           rk = 1x3 numpy array position of other dipole
           k = omega / c
        Returns:
            Ajk element of given position. Should be a 3x3 numpy array
    """
    rjk = np.abs(np.linalg.norm(rj) - np.linalg.norm(rk)) # get rjk
    # get the unit vector rujk
    rujk = np.subtract(rk,rk)/rjk
    # Bracket term of Ajk
    bracketA = k**2*(np.subtract(np.dot(rujk,rujk), np.identity(3))) + (1j*k*rjk - 1)/(rj*k**2)*(np.subtract(3*np.dot(rujk,rujk),np.identity(3)))
    nonBracketA = np.exp(1j*k*rjk)/rjk
    return bracketA + nonBracketA
def complexConjugateGradient(ajk, Einc, N)
    """Returns Polarized Matrix Pk
        Args:
            ajk = 3x3 numpy array from get Ajk
            Einc = Electric field incoming at the j point
            N = number of dipoles
        Normally also would take a guess value but per Draine and Flateu (1994) guess is 0
        Returns:
            Polarized Matrix P numpy array and goal of the project
    """
   # 
