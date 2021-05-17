import numpy as np
import math
import matplotlib.pyplot as plt
import sklearn
from scipy import integrate


####################################################################################################
####################################################################################################
#### Gestion des matrices

def mat_renor(M):
    N=np.copy(M)
    return N/N.max()

def matrice_nan(matrice_distance):
    """
    transforme les 'nan' en 0
    prends en entrée la matrice des distance inter-patch
    """
    M=np.copy(matrice_distance)
    for i in range(0,len(matrice_distance)):
        M[i][i] = 0
    return M

def matrice_stochastique(M):
    """
    Renvoie la version stochastique d'une matrice à coefficients compris entre 0 et 1
    """
    N=np.copy(M)
    for i in range(0,len(M)):
        N[i] = M[i]/np.sum(M[i])
    return N 
    
    
#Finot

def carre(aire,distance):
    longueur=np.sqrt(aire)
    return [[distance-longueur/2,distance+longueur/2],[distance-longueur/2,distance+longueur/2]]


#Kernels

def euclide(b1,b2):
    return np.sqrt((b1[0]-b2[0])**2+(b1[1]-b2[1])**2)
def kernel_2Dt(d,a=100,b=1):
    return b/(np.pi*a)*(1+d**2/a)**(-1-b)
def kernel_exppow(d,a,b):
    return b/(2*np.pi*a**2*math.gamma(2/b))*np.exp(-(d/a)**b)
def kernel_expo(d,a,b=1):
    return 1/(np.pi*a**2)*np.exp(-(d/a)**2)

def kernel_weibull(d,a,b):
    return b/(2*np.pi*a**2)*np.exp(-(b*(d-a)**2/(2*a**2*d)))
def kernel_lognormal(d,a,b):
    return 1/((2*np.pi)**1.5*b*d**2)*np.exp(-np.log((d/a)**2)/(2*b**2))


def kernel_expo_bary(bary1,bary2=[0,0],a=1,b=1):
    return kernel_expo(euclide(bary1,bary2),a,b)
def kernel_2Dt_bary(bary1,bary2=[0,0],a=1,b=1):
    return kernel_2Dt(euclide(bary1,bary2),a,b)
def kernel_weibull_bary(b1,b2=[0,0],a=1,b=1):
    return kernel_weibull(euclide(b1,b2),a,b)
def kernel_lognormal_bary(b1,b2=[0,0],a=1,b=1):
    return kernel_lognormal(euclide(b1,b2),a,b)


#graines

def card_graine_disp_par_patch(compo, aire, taux_disp):
    return compo*aire*taux_disp
def card_graine_restant_par_patch(compo, aire, taux_disp):
    return compo*aire - vect_graine_disp_par_patch(compo,aire,taux_disp)
def card_graine_recu_i_par_j(compo_j,aire_i,aire_j,taux_disp,kernel,d_ij,a,b):
    card_disp = card_graine_disp_par_patch(compo_j, aire_j, taux_disp)
    k = lambda x0,x1 : kernel([x0,x1],[0,0],a,b)
    proba_graine_arrive_i = integrate.nquad(k,carre(aire_i,d_ij))[0]
    return card_disp*proba_graine_arrive_i

def creation_matrice_transition(D,liste_aire,kernel,taux_disp,a,b):
    M=np.zeros((len(liste_aire),len(liste_aire)))
    for i in range(0,len(liste_aire)):
        for j in range(0,len(liste_aire)):
            #proba que graine de j arrive en i
            k = lambda x0,x1 : kernel([x0,x1],[0,0],a,b)
            if i !=j :
                M[i][j] = integrate.nquad(k,carre(liste_aire[i],D[i][j]))[0]
    return M

def matrice_echange(Distance,liste_aire,kernel,rho,pi,a,b):
    M=np.zeros((len(liste_aire)+1,len(liste_aire)+1))
    Proba_kernel = creation_matrice_transition(Distance,liste_aire,kernel,rho,a,b)
    for i in range(0,len(M)-1):
        for j in range(0,len(M)):
            if i == j:
                M[i][j] = (1- rho)*liste_aire[i]
            elif j<(len(M)-1):
                M[i][j] = rho * Proba_kernel[i][j] * liste_aire[j]
            else:
                M[i][j] = pi * liste_aire[i]
    M[-1][-1]=1
    return M