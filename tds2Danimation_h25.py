#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

from vpython import *
import os
import time
import json
import numpy as np
import math
import matplotlib.pyplot as plt


#------------------------------------------------------------------FONCTION_AJOUTER-------------------------------------
# Fonction permettant enregistrer les qdm.

def save_quantite_mouvement():
    # Obtenir le répertoire du script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "p_data.json")  # Chemin complet du fichier

    # Convertir les vecteurs VPython en dictionnaires compréhensibles
    p_dict = [{"x": vect.x, "y": vect.y} for vect in p]

    # Sauvegarder au format JSON
    try:
        with open(filepath, 'w') as f:
            json.dump(p_dict, f, indent=4)  # Sauvegarde
        print(f"Données sauvegardées dans le fichier {filepath}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde : {e}")

#------------------------------------------------------------------FIN_FONCTION_AJOUTER---------------------------------

#------------------------------------------------------------------FONCTION_AJOUTER-------------------------------------
# Fonction permettant enregistrer les distance et le temps entre chaque collision.

def save_distance_temps_collision():
    # Obtenir le répertoire du script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "distance_temps_collision_data.json")  # Chemin complet du fichier

    distance_temps_dict = [{"distance": d, "temps": t} for d, t in zip(distance_collision, temps_collision)]

    # Sauvegarder au format JSON
    try:
        with open(filepath, 'w') as f:
            json.dump(distance_temps_dict, f, indent=4)  # Sauvegarde
        print(f"Données sauvegardées dans le fichier {filepath}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde : {e}")

#------------------------------------------------------------------FIN_FONCTION_AJOUTER---------------------------------

#------------------------------------------------------------------FONCTION_AJOUTER-------------------------------------
# Fonction permettant enregistrer les distance (x et y) et le temps entre chaque collision.

def save_distance_x_y_temps_collision():
    # Obtenir le répertoire du script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(script_dir, "distance_x_y_temps_collision_data.json")  # Chemin complet du fichier

    distance_x_y_temps_dict = [{"distance_x": d_x, "distance_y": d_y, "temps": t} for d_x, d_y, t in zip(distance_collision_x, distance_collision_y, temps_collision)]

    # Sauvegarder au format JSON
    try:
        with open(filepath, 'w') as f:
            json.dump(distance_x_y_temps_dict, f, indent=4)  # Sauvegarde
        print(f"Données sauvegardées dans le fichier {filepath}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde : {e}")

#------------------------------------------------------------------FIN_FONCTION_AJOUTER---------------------------------

#-------------------------------------------------FONCTION_AJOUTER------------------------------------------------------
# Choix de la particule suivit entre chaque collision.

id_particule = 1

#-----------------------------------------------FIN_FONCTION_AJOUTER----------------------------------------------------


# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 200  # change this to have more or fewer atoms
dt = 1E-5  # pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
DIM = 2 #Nombre de degrés de liberté de la simulation 
mass = 4E-3/6E23 # helium mass
Ratom = 0.01 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas( width=750, height=500) # , align='left')
animation.range = L
# animation.title = 'Théorie cinétique des gaz parfaits'
# s = """  Simulation de particules modélisées en sphères dures pour représenter leur trajectoire ballistique avec collisions. Une sphère est colorée et grossie seulement pour l’effet visuel permettant de suivre sa trajectoire plus facilement dans l'animation, sa cinétique est identique à toutes les autres particules.

# """
# animation.caption = s

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d,-d,0), vector(d,-d,0), vector(d,d,0), vector(-d,d,0), vector(-d,-d,0)])

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Atoms = [] # Objet qui contiendra les sphères pour l'animation
p = [] # quantité de mouvement des sphères
apos = [] # position des sphères
pavg = sqrt(2*mass*(DIM/2)*k*T) # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for i in range(Natoms):
    x = L*random()-L/2 # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L*random()-L/2
    z = 0
    if i == 0:  # garde une sphère plus grosse et colorée parmis toutes les grises
        Atoms.append(simple_sphere(pos=vector(x,y,z), radius=0.03, color=color.magenta)) #, make_trail=True, retain=100, trail_radius=0.3*Ratom))\

# -------------------------------------------------FONCTION_AJOUTER------------------------------------------------------
# Met en rouge la particule suivit.

    if i == id_particule:
        Atoms.append(simple_sphere(pos=vector(x,y,z), radius=Ratom, color=color.red))

# -----------------------------------------------FIN_FONCTION_AJOUTER----------------------------------------------------

    else: Atoms.append(simple_sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z)) # liste de la position initiale de toutes les sphères
#    theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # qte de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    p.append(vector(px,py,pz)) # liste de la quantité de mouvement initiale de toutes les sphères

#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions():
    hitlist = []   # initialisation
    r2 = 2*Ratom   # distance critique où les 2 sphères entre en contact à la limite de leur rayon
    r2 *= r2   # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj   # la boucle dans une boucle itère pour calculer cette distance vectorielle dr entre chaque paire de sphère
            if mag2(dr) < r2:   # test de collision où mag2(dr) qui retourne la norme élevée au carré de la distance intersphère dr
                hitlist.append([i,j]) # liste numérotant toutes les paires de sphères en collision
    return hitlist


#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre

#-------------------------------------------------FONCTION_AJOUTER------------------------------------------------------
# Régit le nombre de steps. Init les liste de distance et de temps de collision pour la particule suivit.

nombre_steps = 1000
frame_counter = 0

distance_collision = []
distance_collision_x = []
distance_collision_y = []
temps_collision = []
temps_collision_x_y = []

#-----------------------------------------------FIN_FONCTION_AJOUTER----------------------------------------------------

# -------------------------------------------------------FONCTION_AJOUTER------------------------------------------------
# Fonction permettant de suivre une particule. Determiner la distance et le temps entre chaque collision.

def suivre_particule(id, hitlist, vitesse, step, dt):

    # Initialisation
    position_particule = []
    distance_collision_particule = []
    temps_collision_particule = []
    current_time = step * dt  # Calcule le temps accumulé jusqu'à l'étape

    # Ajoute la position actuelle de la particule suivie
    position_particule.append(apos[id])  # Enregistre la position actuelle de la particule suivie

    # Vérifie si la particule suivie est impliquée dans une collision
    if hitlist:  # S'il y a des collisions détectées
        for paire in hitlist:  # Parcourt les collisions
            i, j = paire  # Indices des particules en collision

            if i == id or j == id:  # Si la particule suivie est impliquée dans une collision
                distance_parcourue = vitesse[id].mag * current_time  # Calcule la distance parcourue
                distance_collision_particule.append(distance_parcourue)  # Stocke la distance parcourue
                temps_collision_particule.append(current_time)  # Stocke le temps écoulé

                current_time = 0  # Réinitialise le temps pour cet intervalle

    return distance_collision_particule, temps_collision_particule


# -----------------------------------------------FIN_FONCTION_AJOUTER---------------------------------------------------

#-------------------------------------------------------FONCTION_AJOUTER------------------------------------------------
# Permet de donner les deplacements x et y. Utile a la question 6 de la partie 1.

def suivre_particule_x_y(id, hitlist, vitesse, step, dt):
    # Initialisation
    position_particule = []
    distance_collision_particule_x = []
    distance_collision_particule_y = []
    temps_collision_particule = []
    current_time = step * dt  # Calcule le temps accumulé jusqu'à l'étape

    # Ajoute la position actuelle de la particule suivie
    position_particule.append(apos[id])  # Enregistre la position actuelle de la particule suivie

    # Vérifie si la particule suivie est impliquée dans une collision
    if hitlist:  # S'il y a des collisions détectées
        for paire in hitlist:  # Parcourt les collisions
            i, j = paire  # Indices des particules en collision

            if i == id or j == id:  # Si la particule suivie est impliquée dans une collision
                distance_parcourue = vitesse[id].mag * current_time  # Calcule la distance parcourue totale

                # Décompose en distance x et distance y
                distance_x = distance_parcourue * (vitesse[id].x / vitesse[id].mag)
                distance_y = distance_parcourue * (vitesse[id].y / vitesse[id].mag)

                distance_collision_particule_x.append(distance_x)  # Stocke la distance x
                distance_collision_particule_y.append(distance_y)  # Stocke la distance y
                temps_collision_particule.append(current_time)  # Stocke le temps écoulé

                current_time = 0  # Réinitialise le temps pour cet intervalle

    return distance_collision_particule_x, distance_collision_particule_y, temps_collision_particule

#---------------------------------------------------FIN_FONCTION_AJOUTER------------------------------------------------


#-------------------------------------------------------FONCTION_AJOUTER------------------------------------------------
# Remplacement de la While par une boucle for qui prend en compte un nombre de Steps.

for step in range(nombre_steps):
    current_time = step * dt  # Update current time based on the step

    rate(300)  # Control the simulation speed

        # Save quantities every 1000 iterations
    if frame_counter % 1000 == 0:
        save_quantite_mouvement()

    frame_counter += 1

#---------------------------------------------------FIN_FONCTION_AJOUTER------------------------------------------------

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []   # vitesse instantanée de chaque sphère
    deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for i in range(Natoms):
        vitesse.append(p[i]/mass)   # par définition de la quantité de nouvement pour chaque sphère
        deltax.append(vitesse[i] * dt)   # différence avant pour calculer l'incrément de position
        Atoms[i].pos = apos[i] = apos[i] + deltax[i]  # nouvelle position de l'atome après l'incrément de temps dt

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE ####
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2:
            if loc.x < 0: p[i].x =  abs(p[i].x)  # renverse composante x au mur de gauche
            else: p[i].x =  -abs(p[i].x)   # renverse composante x au mur de droite
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)  # renverse composante y au mur du bas
            else: p[i].y =  -abs(p[i].y)  # renverse composante y au mur du haut


    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist = checkCollisions()

# -------------------------------------------------FONCTION_AJOUTER-----------------------------------------------------
    # Permet de détecter les collisions et de créer des listes avec la distance parcourue et le temps entre les collision.

    distance_collision_particule, temps_collision_particule = suivre_particule(id_particule, hitlist, vitesse, step, dt)

    # Ajouter les données collectées
    distance_collision.extend(distance_collision_particule)
    temps_collision.extend(temps_collision_particule)

    # Enregistre les donnees de distance.
    save_distance_temps_collision()

    # Permet de détecter les collisions et de créer des listes avec la distance (x et y) parcourue et le temps entre les collision.
    distance_x, distance_y, temps_collision_x_y= suivre_particule_x_y(id_particule, hitlist, vitesse, step, dt)

    # Ajouter les donnees collectees
    distance_collision_x.extend(distance_x)
    distance_collision_y.extend(distance_y)
    temps_collision_x_y.extend(temps_collision_x_y)

    save_distance_x_y_temps_collision()


# -----------------------------------------------FIN_FONCTION_AJOUTER---------------------------------------------------

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for ij in hitlist:

        # définition de nouvelles variables pour chaque paire de sphères en collision
        i = ij[0]  # extraction du numéro des 2 sphères impliquées à cette itération
        j = ij[1]
        ptot = p[i]+p[j]   # quantité de mouvement totale des 2 sphères
        mtot = 2*mass    # masse totale des 2 sphères
        Vcom = ptot/mtot   # vitesse du référentiel barycentrique/center-of-momentum (com) frame
        posi = apos[i]   # position de chacune des 2 sphères
        posj = apos[j]
        vi = p[i]/mass   # vitesse de chacune des 2 sphères
        vj = p[j]/mass
        rrel = posi-posj  # vecteur pour la distance entre les centres des 2 sphères
        vrel = vj-vi   # vecteur pour la différence de vitesse entre les 2 sphères

        # exclusion de cas où il n'y a pas de changements à faire
        if vrel.mag2 == 0: continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        if rrel.mag > Ratom: continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta) où theta is the angle between vrel and rrel:
        dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        alpha = asin(dy/(2*Ratom))  # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
        deltat = d/vrel.mag         # time spent moving from first contact to position inside atom

        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        posi = posi-vi*deltat   # back up to contact configuration
        posj = posj-vj*deltat
        pcomi = p[i]-mass*Vcom  # transform momenta to center-of-momentum (com) frame
        pcomj = p[j]-mass*Vcom
        rrel = hat(rrel)    # vecteur unitaire aligné avec rrel
        pcomi = pcomi-2*dot(pcomi,rrel)*rrel # bounce in center-of-momentum (com) frame
        pcomj = pcomj-2*dot(pcomj,rrel)*rrel
        p[i] = pcomi+mass*Vcom # transform momenta back to lab frame
        p[j] = pcomj+mass*Vcom
        apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération
        apos[j] = posj+(p[j]/mass)*deltat


