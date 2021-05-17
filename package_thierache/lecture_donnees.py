import pandas as pd
import numpy as np

def liste_PA(PA_sheet):
    PA_mat = PA_sheet.values
    PA=[]
    for i in range(0,len(PA_mat[0,:])):
        PA.append(sum(PA_mat[:,i]))
    return PA

def extraction_S(PA_value):
    S=0
    for i in range(0,len(PA_value)):
        if sum(PA_value[i,1:])>0:
            S+=1
    return S

def creation_liste_disp(nom_espece_total, nom_espece_local, disp_total):
    """
    crée la liste des tupe de dispertion LOCALE
    """
    disp_local = []
    for i in range(0,len(nom_espece_local)):
        for j in range(0,len(nom_espece_total)):
            if nom_espece_local[i]==nom_espece_total[j]:
                disp_local.append(disp_total[j])
    return disp_local

def barycentre_a_matrice(liste_bary):
    M=np.zeros((len(liste_bary),len(liste_bary)))
    for i in range(0,len(liste_bary)):
        for j in range(0,len(liste_bary)):
            M[i][j]=np.sqrt((liste_bary[i][0]-liste_bary[j][0])**2+(liste_bary[i][1]-liste_bary[j][1])**2)
    return M

def extraction_data_OB(df_PA, df_aire, df_disp, df_barycentre, OB=True, mat=0):
    #Fichier Présence Absence, Aire, Dispersion ET BARYCENTRE
    PA=pd.read_excel(df_PA)
    aire = pd.read_excel(df_aire)
    disp = pd.read_excel(df_disp)
    bary = pd.read_excel(df_barycentre)
    
    if OB == True:
        liste_xls = [0,2,4]
        S = len(pd.read_excel(df_PA, sheet_name=0, index_col=0)['OB_01'].tolist())
    else:
        liste_xls = [1,3,5]
        S = len(pd.read_excel(df_PA, sheet_name=1, index_col=0)['OT_02'].tolist())

    type_dispersion = pd.read_excel(df_disp, sheet_name=0, index_col=0)
    
    PA_OB = pd.read_excel(df_PA, sheet_name=liste_xls[0], index_col=0)

    liste_aire = []
    for i in liste_xls:
        liste_aire.append(pd.read_excel(df_aire, sheet_name=i, index_col=0))
    
    #Du barycentre à la matrice
    liste_matrice = []
    liste_bary=[]
    for i in liste_xls:
        coord_bary = pd.read_excel(df_barycentre, sheet_name=i, index_col=0)
        liste_bary.append(coord_bary.values[:,0:])
        if mat==0:
            liste_matrice.append( barycentre_a_matrice(coord_bary.values[:,0:]))
        else:
            matrice = pd.read_excel(mat)
            liste_matrice.append(pd.read_excel(mat, sheet_name=i, index_col=0))
    
    
    #Par époque : [aire, matrice]
    OB_ajd = [liste_aire[0], liste_matrice[0],liste_bary[0]]
    OB_5060 = [liste_aire[1], liste_matrice[1],liste_bary[1]]
    OB_EM = [liste_aire[2], liste_matrice[2],liste_bary[2]]
    #R_reel
    R_reel = liste_PA(PA_OB)
    
    #l_disp_local
    l_disp = type_dispersion
    #aire_ajd, matrice_ajd = OB_ajd
    #liste_aire_ajd = aire_ajd['Shape_Area_m2'].tolist()
    
    liste_espece_beauvaisis = PA['Species names following BDNFF v 2010.01.05 Tela Botanica'].tolist()
    l_espece_local = liste_espece_beauvaisis
    type_disp = l_disp['Unnamed: 1'].tolist()
    nom_espece_total=l_disp.index
    l_disp_local =  creation_liste_disp(nom_espece_total, liste_espece_beauvaisis,type_disp)
    
    
    return S, R_reel, OB_EM, OB_5060, OB_ajd, l_espece_local, l_disp_local, PA
    
def extraction_data(df_PA='PA.xlsx',df_disp='disp.xlsx', landscape='OC', mat=0):
    #Fichier Présence Absence, Aire, Dispersion ET BARYCENTRE
    PA=pd.read_excel(df_PA)
    disp = pd.read_excel(df_disp)
    df_aur = pd.read_excel('DataAurelien.xls')
    
    if landscape == 'OC':
        i_min, i_max = 184, 214
        S = extraction_S(pd.read_excel(df_PA, sheet_name=2).values)
        
    type_dispersion = pd.read_excel(df_disp, sheet_name=0, index_col=0)
    
    PA_OB = pd.read_excel(df_PA, sheet_name=2, index_col=0)
    
    OX = df_aur.iloc[i_min:i_max,:]
    liste_aire = OX['Surface'].values
    
    x_bary, y_bary = OX['CentroidX'].values, OX['CentroidY'].values
    
    #Du barycentre à la matrice
    liste_bary=[]
    for i in range(0,len(x_bary)):
        liste_bary.append([x_bary[i],y_bary[i]])
        
    matrice=barycentre_a_matrice(liste_bary)
        
    #Par époque : [aire, matrice]
    OB_ajd = [liste_aire, matrice,liste_bary]
    
    #R_reel
    R_reel = liste_PA(PA_OB)
    
    #l_disp_local
    l_disp = type_dispersion
    
    liste_espece_beauvaisis = PA['Species names following BDNFF v 2010.01.05 Tela Botanica'].tolist()
    l_espece_local = liste_espece_beauvaisis
    type_disp = l_disp['Unnamed: 1'].tolist()
    nom_espece_total=l_disp.index
    l_disp_local =  creation_liste_disp(nom_espece_total, liste_espece_beauvaisis,type_disp)
    
    
    return S, R_reel, OB_ajd, l_espece_local, l_disp_local
    
    
    
def extraction_fragment_virtuel(df_virtuel):
    virtuel = pd.read_excel(df_virtuel, sheet_name = 1, index_col=0)
    aire = virtuel['Shape_Area_m2'].tolist()
    richesse = virtuel['nb_especes']
    
    return aire, richesse
    
#def clean_patch