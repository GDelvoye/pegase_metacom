import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, "/home/delvoye/Documents/Jupyter/Projets/chap3/package_thierache")
import initialisation as ini
import matrice as mat
import analyse_resultat as ana
import sauvegarde as svg
import pandas as pd

def liste_Gene(l_disp):
    l=[]
    for i in range(0,len(l_disp)):
        if l_disp[i]=='G':
        #if l_disp[i]!='S':
            l.append(i)
    return l
    
def liste_Spe(l_disp):
    l=[]
    for i in range(0,len(l_disp)):
        if l_disp[i]=='S':
            l.append(i)
    return l

class Espece():
    def __init__(self,nb_graine=1,disp_graine=1,nom=''):
        self.nb_graine=nb_graine
        self.disp_graine = disp_graine  
        self.nom = nom
        

class Patch():
    def __init__(self,aire,nom,barycentre=np.array([0,0])):
        self.aire = aire
        #self.composition = composition  
        self.nom = nom
        self.barycentre = barycentre
        

class Metacom():
    def __init__(self,liste_patch,S, tableau_composition, nom_especes,label, matrice_distance, liste_aire,\
                 pi=5e-6,rho=0.1,kernel=mat.kernel_expo_bary,fenetre = 'OB',a=10000,b=1,aS=1,bS=1,rhoS=0.01,\
                 type_disp=None,nbS=1,s=10,delta_s=5, pi_s = 6e-06):
        self.liste_patch = liste_patch
        self.S = S
        self.type_disp = label
        self.composition = tableau_composition
        self.nom_especes = nom_especes
        self.liste_aire = list(liste_aire)
        self.matrice_distance=matrice_distance
        self.histoire_composition=[]
        self.pool = np.ones(S)/S
        self.pi=pi
        self.pi_s = pi_s
        self.s = s
        self.delta_s=delta_s
        self.rho = rho
        self.rhoS = rhoS
        self.kernel = kernel
        self.a=a
        self.aS=aS
        self.fenetre = fenetre
        self.b=b
        self.bS=bS
        self.nbS=nbS
        self.MS= mat.matrice_echange(self.matrice_distance,\
                                     self.liste_aire,self.kernel,self.rhoS,\
                                     self.pi_s,self.aS,self.bS)
        self.M = mat.matrice_echange(self.matrice_distance,\
                                     self.liste_aire,self.kernel,self.rho,self.pi,\
                                     self.a,self.b)
        self.ind_G = liste_Gene(self.type_disp)
        self.ind_S = liste_Spe(self.type_disp)
        
    def richesse_observee(self):
        if self.fenetre == 'OB':
            df = pd.read_excel('PA_OB.xls')
        elif self.fenetre == 'OT':
            df = pd.read_excel('PA_OT.xls')
        else:
            df = pd.read_excel('PA_OC.xls')
        R_obs = []
        for i in range(3,len(df.values[0])):
            R_obs.append(sum(df.values[:,i]))
        return R_obs
    
    def stats(self, full = False):
        return ana.stats_metacom(self.composition[:-1,], full)
    def plot_territoire(self):
        ana.plot_territoire(self)
    def plot_territoire2(self):
        ana.plot_territoire2(self)
    def plot_diversite(self, show = 'obs'):
        ana.plot_diversite(self, show = show)
    
    
    