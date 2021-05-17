import pickle

def recuperation(fichier):
    with open(fichier, 'rb') as nom_fichier:
        depickler = pickle.Unpickler(nom_fichier)
        metacom=depickler.load()
    return metacom

def sauvegarde(metacom, nom_fichier):
    with open(nom_fichier, 'wb') as fichier:
        pickler = pickle.Pickler(fichier)
        pickler.dump(metacom)