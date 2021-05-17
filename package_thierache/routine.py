import numpy as np
import random as rd
import sys
sys.path.insert(0, "/home/delvoye/Documents/Jupyter/Projets/chap3/package_thierache")

import matrice as mat
import classes as cl
import wf_discret as wf
import analyse_resultat as ana
import initialisation as ini
import lecture_donnees as lecture



def histo_compo(metacom):
    liste_compo=[]
    for i in range(0,len(metacom.liste_patch)):
        liste_compo.append(metacom.liste_patch[i].composition)
    metacom.histoire_compo.append(liste_compo)


def ajout_patch(metacom, nouveau_bary, a, compo=False):
    S = metacom.S
    c=np.zeros(S)
    c[0]=1
    nouveau_patch = cl.Patch(a,"nouveau patch",nouveau_bary)
    metacom.liste_patch.append(nouveau_patch)
    metacom.liste_aire.append(np.array([a]))   
    

def augmentation_S(tableau_compo,ajout_S=100):
    tab=np.copy(tableau_compo)
    tab=np.concatenate((tableau_compo,np.zeros((len(tab),ajout_S))),axis=1)
    tab[-1]=np.ones(len(tableau_compo[0])+ajout_S)/(len(tableau_compo[0])+ajout_S)
    return tab

def modif_composition(m, liste_patch):
    compo_nouveau = np.zeros(len(m.composition[0]))
    compo_nouveau[0]=1 #fixation de l'espèce 1 dans les nouveaux patch
    pool = m.composition[-1]
    compo = list(m.composition[:-1])
    for i in range(0,len(liste_patch)):
        compo.append(compo_nouveau)
    compo.append(pool)
    m.composition = np.array(compo)

def modif_metacom(metacom, a=None,b=None, pi=None, a_S =None, ker=None, new_patch=None,\
                  modif_l_aire = None, rho=None, aS=None):
    #I. Vérification des attributs
    if modif_l_aire !=None:
        for i in range(0,len(modif_l_a)):
            metacom.liste_aire[modif_l_aire[i][0]] = modif_l_aire[i][1]
            metacom.liste_patch[modif_l_aire[i][0]].aire = modif_l_aire[i][1]
    if a!= None:
        metacom.a = a
    if b != None:
        metacom.b = b
    if pi!=None:
        metacom.pi = pi
    if a_S != None:
        metacom.S+=a_S
        metacom.pool = np.ones(metacom.S)/sum(np.ones(metacom.S))
        metacom.composition=augmentation_S(metacom.composition,a_S)
        for i in range(0,a_S):
            if rd.random()< 175./245:
                metacom.type_disp.append('G')
            else:
                metacom.type_disp.append('S')
    
    if ker != None:
        metacom.kernel = ker
    if rho!= None:
        metacom.rho = rho       

    
    l_bary=[]
    for patch in metacom.liste_patch:
        l_bary.append(patch.barycentre)
    if new_patch!=None:
        for patch in new_patch:
            ajout_patch(metacom,patch[0],patch[1])
            l_bary.append(patch[0])
            modif_composition(metacom,new_patch)
    metacom.matrice_distance=lecture.barycentre_a_matrice(l_bary)
    metacom.M=mat.matrice_echange(metacom.matrice_distance,metacom.liste_aire,\
                                                 metacom.kernel,metacom.rho, \
                                  metacom.pi,metacom.a,metacom.b)
    
    

def routine_meta(metacom, pas=100,indiv_par_m2=1,\
                    pool=True, historique = False): 
    for i in range(0,pas):
        wf.pas_neutre(metacom, indiv_par_m2,pool)
        if historique == True:
            metacom.histoire_composition.append(metacom.composition)
    if historique == False:
        metacom.histoire_composition.append(metacom.composition)



def routine_present(liste_nom_fichier,pas=200,s=-1,theorique=True,indiv_par_m2=1, neutre=True,\
                    OB=True,epoque='ajd',pool=True,tx_pool=1e-5,\
                    ker = mat.kernel_2Dt_bary,par_a=1,par_b=1,tx_dispersion = 0.1,lisiere=0):
    
    #I.Extraction fichier
    S, R_reel, OB_EM, OB_5060, OB_ajd, l_espece_local, l_disp_local, PA = lecture.extraction_data_OB(liste_nom_fichier[0],liste_nom_fichier[1],liste_nom_fichier[2],liste_nom_fichier[3],OB,mat=lisiere)
    #II. Choix de l'époque
    if epoque =='ajd':
        aire, matrice, bary = OB_ajd
    elif epoque == '5060':
        aire, matrice, bary = OB_5060
    else:
        aire, matrice, bary = OB_EM
        
    liste_aire = aire.values
    if lisiere!=0:
        matrice = mat.matrice_nan(matrice.values)
    
    #III. Initialisation de liste_patch
    if type(s) == list:
        #On donne une liste de patchs en entrée
        liste_patch = []
        for i in range(0,len(liste_aire)):
            liste_patch.append(cl.Patch(int(liste_aire[i]),s[i],aire.index[i],bary[i]))
    
    elif s == - 1:
        liste_patch = []
        for i in range(0,len(liste_aire)):
            if theorique==True:
                #On ini avec la richesse théorique
                s_theo = ini.richesse_theorique(liste_aire[i])
                liste_patch.append(cl.Patch(int(liste_aire[i]),\
                                            ini.mix_partiel(S,s_theo),aire.index[i],bary[i]))
            else:
                #On ini avec les PA observées
                liste_patch.append(cl.Patch(int(liste_aire[i]),\
                                            np.array(list(PA.values[:,i+1]/sum(PA.values[:,i+1]))),aire.index[i],\
                                           bary[i]))
    elif s=="reel":
        liste_patch = []
        for i in range(0,len(liste_aire)):
            liste_patch.append(cl.Patch(int(liste_aire[i]),ini.mix_partiel(S,R_reel[i]),aire.index[i],bary[i]))
    
    else:
        #On ini avec s espèces par patchs en équirépartition
        liste_patch = []
        for i in range(0,len(liste_aire)):
            liste_patch.append(cl.Patch(int(liste_aire[i]),ini.mix_partiel(S,s),aire.index[i],bary[i]))

    #IV. Création de la métacommunauté
    metacom = cl.Metacom(liste_patch, S, l_espece_local, l_disp_local,\
                         matrice,liste_aire,taux_pool=tx_pool,\
                         kernel=ker,taux_disp=tx_dispersion,a=par_a,b=par_b)
    histo_compo(metacom)
    metacom.histoire_R.append(ana.richesse_specifique_liste(metacom.histoire_compo[-1]))
    metacom.histoire_gamma.append(ana.diversite_gamma(metacom.histoire_compo[-1]))
    
    #V. Simulation forward
    for i in range(0,pas):
        if neutre == False:
            wf.pas(metacom, indiv_par_m2,pool)
        else:
            wf.pas_neutre(metacom, indiv_par_m2,pool)
        histo_compo(metacom)
        metacom.histoire_R.append(ana.richesse_specifique_liste(metacom.histoire_compo[-1]))
        metacom.histoire_gamma.append(metacom.histoire_compo[-1])
        
    return metacom




