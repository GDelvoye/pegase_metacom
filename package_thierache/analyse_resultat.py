import numpy as np
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
import pandas as pd
from matplotlib import colors
import sauvegarde as svg
import routine as rt




#Diversité

def richesse_commu(compo,ind_disp=[]):
    R=0
    for i in range(0,len(compo)):
        if ind_disp == []:
            if compo[i]>0:
                R+=1
        else:
            if (i in ind_disp) and compo[i]>0:
                R+=1
    #for j in compo:
     #   if j >0.:
      #      R+=1
    return R

##### beta
def ana_beta_pair(c1,c2):
    R1 = richesse_commu(c1)
    R2 = richesse_commu(c2)
    a,b,c = 0, 0, 0#commun, unique pauvre, unique riche
    for i in range(0,len(c1)):
        if c1[i]>0:
            if c2[i]>0:
                a+=1
            else: #1 oui, 2 non
                if R1>R2:#1 est le riche
                    c+=1
                else:
                    b+=1
        elif c2[i]>0:#2 oui, 1 non
            if R1>R2:#1 est le riche
                b+=1
            else:
                c+=1
    return a, b, c

def ana_beta_pairVRAI(c1,c2):
    R1 = richesse_commu(c1)
    R2 = richesse_commu(c2)
    a,b,c = 0, 0, 0#commun, unique pauvre, unique riche
    for i in range(0,len(c1)):
        if c1[i]>0:
            if c2[i]>0:
                a+=1
            else: #1 oui, 2 non
                b+=1
        elif c2[i]>0:#2 oui, 1 non
                c+=1
    return a, b, c

def beta_pair(c1,c2,typ="sorensen"):
    a, b, c = ana_beta_pair(c1,c2)
    if typ == "sorensen":
        beta_nested = (c-b)/(2*a+b+c)*(a/(a+b))
        beta_turn = b/(b+a)
        return beta_nested+beta_turn, beta_turn, beta_nested
    else:
        beta_nested = (c-b)/(a+b+c)*(a/(a+2*b))
        beta_turn = 2*b/(2*b+a)
        return beta_nested+beta_turn, beta_turn, beta_nested


def diversite_gamma(liste_compo):
    R=0
    for s in range(0,len(liste_compo[0])):
        R_s = 0
        for i in range(0,len(liste_compo)):
            R_s+=liste_compo[i][s]
        if R_s >0:
            R+=1
    return R

def beta_multi(liste_compo,typ="sorensen"):
    S_T = diversite_gamma(liste_compo)
    Sum_Si = 0
    Sum_min, Sum_max = 0,0
    for i in range(0,len(liste_compo)):
        Sum_Si += richesse_commu(liste_compo[i])
        for j in range(0,len(liste_compo)):
            if i < j:
                a,b,c = ana_beta_pair(liste_compo[i],liste_compo[j])
                Sum_min+=b
                Sum_max+=c
    BETA_SIM = Sum_min/((Sum_Si - S_T)+Sum_min)
    BETA_SNE = (Sum_max - Sum_min)/(2*(Sum_Si - S_T) + Sum_min + Sum_max) * (Sum_Si - S_T)/(Sum_Si - S_T + Sum_min)
    return (BETA_SIM + BETA_SNE,BETA_SIM, BETA_SNE)
def beta_multi(liste_compo,typ="sorensen"):
    S_T = diversite_gamma(liste_compo)
    Sum_Si = 0
    Sum_min, Sum_max = 0,0
    for i in range(0,len(liste_compo)):
        Sum_Si += richesse_commu(liste_compo[i])
        for j in range(0,len(liste_compo)):
            if i < j:
                a,b,c = ana_beta_pair(liste_compo[i],liste_compo[j])
                Sum_min+=min(b,c)
                Sum_max+=max(b,c)
    BETA_SIM = Sum_min/((Sum_Si - S_T)+Sum_min)
    BETA_SNE = (Sum_max - Sum_min)/(2*(Sum_Si - S_T) + Sum_min + Sum_max) * (Sum_Si - S_T)/(Sum_Si - S_T + Sum_min)
    return (BETA_SIM + BETA_SNE,BETA_SIM, BETA_SNE)
                

def richesse_specifique_liste(liste_commu, ind_disp =[]):
    #Retourne le R d'une liste de patchs
    l_R = []
    for i in range(0,len(liste_commu)):
        l_R.append(richesse_commu(liste_commu[i],ind_disp))
    return l_R

def stats_metacom(liste_commu, full=False):
    #renvoie alpha, beta, gamma
    if full == True:
        return richesse_specifique_liste(liste_commu), beta_multi(liste_commu), diversite_gamma(liste_commu)
    else:
        return richesse_specifique_liste(liste_commu), beta_multi(liste_commu)[0], diversite_gamma(liste_commu)

def plot_territoire(metacom,size = (8,6),area_expo = 1.,area_deno=500,vmin=None,vmax=None,alpha=1,cmap='Greens'):
    l_x = []
    l_y = []
    l_s = []
    l_c = []
    r_obs=metacom.richesse_observee()
    a,b,g = stats_metacom(metacom.composition[:-1,])
    for i in range(0,len(metacom.liste_patch)):
        l_x.append(metacom.liste_patch[i].barycentre[0])
        l_y.append(metacom.liste_patch[i].barycentre[1])
        l_s.append((metacom.liste_patch[i].aire/area_deno)**area_expo)
        if len(r_obs)==len(metacom.liste_patch):
            l_c.append(a[i]-r_obs[i])
    if vmin == None:
        if len(r_obs)==len(metacom.liste_patch):
            vmin = min(l_c)
            vmax=max(l_c)
        else:
            vmin = 0
            vmax = 1
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    plt.figure(figsize=size)
    sc=plt.scatter(l_x,l_y,s=l_s,c=l_c,marker='o',cmap=cmap,norm=norm,alpha=alpha,linewidths=0.7,edgecolors='black')
    plt.colorbar(sc)
    #plt.xlim(0,5000)
    #plt.ylim(0,5000)
    plt.show()
    
def plot_territoire2(metacom,size = (8,8),area_expo = 1.,area_deno=500,vmin=None,vmax=None,alpha=1,cmap='Greens'):
    l_x = []
    l_y = []
    l_s = []
    l_c = []
    r_obs=metacom.richesse_observee()
    a,b,g = stats_metacom(metacom.composition[:-1,])
    for i in range(0,len(metacom.liste_patch)):
        l_x.append(metacom.liste_patch[i].barycentre[0])
        l_y.append(metacom.liste_patch[i].barycentre[1])
        l_s.append((metacom.liste_patch[i].aire/area_deno)**area_expo)
        if i >= 30:
            l_c.append(0)
        else:
            l_c.append(1)
    if vmin == None:
        if len(r_obs)==len(metacom.liste_patch):
            vmin = min(l_c)
            vmax=max(l_c)
        else:
            vmin = 0
            vmax = 1
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    plt.figure(figsize=size)
    sc=plt.scatter(l_x,l_y,s=l_s,c=l_c,marker='o',norm=norm,alpha=alpha,linewidths=0.7,edgecolors='black')
    
    #plt.xlim(0,5000)
    #plt.ylim(0,5000)
    plt.show()

def plot_diversite(metacom,size=6, show='full'):
    if metacom.fenetre == 'OB':
        #nb_patch_ori = 30
        gamma_obs = 192
    elif metacom.fenetre =='OC':
        gamma_obs = 204
        #nb_patch_ori = 30
    else:#OT
        gamma_obs = 196
        #nb_patch_ori = 29
    R_obs = metacom.richesse_observee()
    R_sim = richesse_specifique_liste(metacom.composition[:-1,])
    R_sim_G = richesse_specifique_liste(metacom.composition[:-1,], metacom.ind_G)
    R_sim_S = richesse_specifique_liste(metacom.composition[:-1,], metacom.ind_S)
    
    plt.figure(figsize=(size,size))
    plt.plot(R_obs,'rx',label='observed')        
    print("Obs VS Sim: "+str(mse(R_obs,R_sim)))
    plt.plot(R_sim,'go',label='simulated')
    if show=='full':
        plt.plot(R_sim_G,'bo',label='generalist')
        plt.plot(R_sim_S,'yo',label='specialist')
    a,b,g=metacom.stats()
    plt.title("Gamma : "+str(g)+"/"+str(gamma_obs)+", a ="+str(metacom.a)+", pi ="+str(metacom.pi)+\
              ", rho ="+str(metacom.rho))
    plt.legend()
    
    plt.show()

    
def monte_carlo_diversite(fichier, n=10, T=10000):
    #renvoie les 3 listes de diversités de taille n
    la, lb, lg =[], [], []
    for i in range(0,n):
        m=svg.recuperation(fichier)
        rt.routine_meta(m, T, pool=True)
        a, b, g=stats_metacom(m.composition[:-1])
        la.append(a)
        lb.append(b)
        lg.append(g)
    return la, lb, lg


  
def mse(l1,l2):
    #MSE
    mse=0
    for i in range(0,len(l1)):
        mse+=(l1[i]-l2[i])**2
    return mse
