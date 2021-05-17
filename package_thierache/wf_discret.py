import numpy as np
import sys
sys.path.insert(0, "/home/delvoye/Documents/Jupyter/Projets/chap3/package_thierache")
import matrice as mat
import initialisation as ini

########################################
def wf_communaute(aire, composition, indiv_par_m2=1):
    return np.random.multinomial(aire/indiv_par_m2, composition)/(aire/indiv_par_m2)

def wf(metacom, indiv_par_m2=1):
    for i in range(0,len(metacom.liste_patch)):
        metacom.composition[i] = wf_communaute(\
                        metacom.liste_patch[i].aire,metacom.composition[i]/sum(metacom.composition[i]), indiv_par_m2)


#Modèle neutre

def echange(metacom, pool = True):
    if pool == True:
        metacom.composition = np.dot(metacom.M, metacom.composition)
    else:
        metacom.composition[:-1] = np.dot(metacom.M[:-1,:-1], metacom.composition[:-1])
    
def pas_neutre(metacom, indiv_par_m2=1,pool=True):
    echange(metacom,pool)
    wf(metacom, indiv_par_m2)
    
#Modèle niche
def echange_spe(metacom, pool = True):
    compo_G = np.dot(metacom.M, metacom.composition[:,metacom.ind_G])
    compo_S = np.dot(metacom.MS, metacom.composition[:,metacom.ind_S]*metacom.nbS)
    metacom.composition[:,metacom.ind_G] = compo_G
    metacom.composition[:,metacom.ind_S] = compo_S

def wf_communaute_sel(aire, composition, indiv_par_m2=1):
    return np.random.multinomial(aire/indiv_par_m2, composition)/(aire/indiv_par_m2)

def selection(communaute,s,ind_S):
    c=[]
    for i in range(0,len(communaute)):
        if i in ind_S:
            c.append((1+s)*communaute[i])
        else:
            c.append(communaute[i])
    return c
        

def wf_sel(metacom, indiv_par_m2=1):
    ind_sel=[0,4,5,8]
    for i in range(0,len(metacom.liste_patch)):
        if i not in ind_sel:
            comm = selection(metacom.composition[i],metacom.s/metacom.delta_s,metacom.ind_S)
            #metacom.composition[i] = wf_communaute(\
                #metacom.liste_patch[i].aire,metacom.composition[i]/sum(metacom.composition[i]), indiv_par_m2)
            metacom.composition[i] = wf_communaute(\
                        metacom.liste_patch[i].aire,comm/sum(comm), indiv_par_m2)
        else:
            comm = selection(metacom.composition[i],metacom.s,metacom.ind_S)
            metacom.composition[i] = wf_communaute(\
                        metacom.liste_patch[i].aire,comm/sum(comm), indiv_par_m2)
        
            
            

    
def pas_spe(metacom, indiv_par_m2=1,pool=True):
    echange_spe(metacom,pool)
    wf_sel(metacom, indiv_par_m2)