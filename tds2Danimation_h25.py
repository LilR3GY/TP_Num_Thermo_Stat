#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt
import json
import os

def save_quantite_mouvement():
    """
    Sauvegarde la quantité de mouvement actuelle (liste `p`) dans un fichier JSON
    dans le même répertoire que ce fichier script.
    """
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

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 2  # change this to have more or fewer atoms
dt = 1e-5  # pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
DIM = 2  # Nombre de degrés de liberté de la simulation
mass = 4e-3 / 6e23  # helium mass
Ratom = 0.01  # wildly exaggerated size of an atom
k = 1.4e-23  # Boltzmann constant
T = 300  # around room temperature

#### CANEVAS DE FOND ####
L = 1  # container is a cube L on a side
gray = color.gray(0.7)  # color of edges of container and spheres below
animation = canvas(width=750, height=500)  # , align='left')
animation.range = L
# animation.title = 'Théorie cinétique des gaz parfaits'
# s = """  Simulation de particules modélisées en sphères dures pour représenter leur trajectoire ballistique avec collisions. Une sphère est colorée et grossie seulement pour l’effet visuel permettant de suivre sa trajectoire plus facilement dans l'animation, sa cinétique est identique à toutes les autres particules.

# """
# animation.caption = s

#### ARÊTES DE BOÎTE 2D ####
d = L / 2 + Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append(
    [
        vector(-d, -d, 0),
        vector(d, -d, 0),
        vector(d, d, 0),
        vector(-d, d, 0),
        vector(-d, -d, 0),
    ]
)

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Atoms = []  # Objet qui contiendra les sphères pour l'animation
p = []  # quantité de mouvement des sphères
apos = []  # position des sphères
pavg = sqrt(2 * mass * (DIM / 2) * k * T)  # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for i in range(Natoms):
    x = (
        L * random() - L / 2
    )  # position aléatoire qui tient compte que l'origine est au centre de la boîte
    y = L * random() - L / 2
    z = 0
    if i == 0:  # garde une sphère plus grosse et colorée parmis toutes les grises
        Atoms.append(
            simple_sphere(pos=vector(x, y, z), radius=0.03, color=color.magenta)
        )  # , make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else:
        Atoms.append(simple_sphere(pos=vector(x, y, z), radius=Ratom, color=gray))
    apos.append(vec(x, y, z))  # liste de la position initiale de toutes les sphères
    #    theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2 * pi * random()  # direction aléatoire pour la quantité de mouvement
    px = pavg * cos(phi)  # qte de mvt initiale selon l'équipartition
    py = pavg * sin(phi)
    pz = 0
    p.append(vector(px, py, pz))  # liste de la quantité de mouvement initiale de toutes les sphères


#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions():
    hitlist = []  # initialisation
    r2 = (2 * Ratom)  # distance critique où les 2 sphères entre en contact à la limite de leur rayon
    r2 *= r2  # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i):
            aj = apos[j]
            dr = (ai - aj)  # la boucle dans une boucle itère pour calculer cette distance vectorielle dr entre chaque paire de sphère
            if (mag2(dr) < r2):  # test de collision où mag2(dr) qui retourne la norme élevée au carré de la distance intersphère dr
                hitlist.append([i, j])  # liste numérotant toutes les paires de sphères en collision
    return hitlist  # retourne la liste des paires de sphères en collision

def follow_the_bouncing_ball():
    """
    Fonction pour suivre la trajectoire d'une sphère particulière dans l'animation
    et calculer le temps et la distance entre chaque collision, en ignorant les
    collisions avec les murs.
    """
    # Crée une sphère à suivre
    sphere_to_follow = simple_sphere(pos=vector(0, 0, 0), radius=0.03, color=color.cyan)
    print(1)
    # Initialisation des variables pour le calcul du temps et de la distance entre les collisions
    previous_collision_time = 0
    collision_times = []
    collision_distances = []

    simulation_time = 0
    max_simulation_time = 1  # Suivre la sphère pendant les 10 premières secondes

    while running and simulation_time < max_simulation_time:
        rate(300)  # Limite la vitesse de calcul de la simulation

        # Met à jour la position de la sphère à suivre
        sphere_to_follow.pos = Atoms[0].pos

        # Vérifie les collisions
        hitlist = checkCollisions()
        print(2)
        # Calcule le temps et la distance entre les collisions pour la sphère à suivre
        for ij in hitlist:
            if 0 in ij:  # Si la sphère à suivre est impliquée dans une collision
                current_time = simulation_time
                collision_time = current_time - previous_collision_time
                previous_collision_time = current_time

                # Calcule la distance parcourue par la sphère à suivre
                i, j = ij
                posi = apos[i]
                posj = apos[j]
                distance = mag(posi - posj)

                # Ajoute les données de collision aux listes
                collision_times.append(collision_time)
                collision_distances.append(distance)

        # Mise à jour de la position des sphères
        for i in range(Natoms):
            vitesse = p[i] / mass
            deltax = vitesse * dt
            Atoms[i].pos = apos[i] = apos[i] + deltax
        print(3)
        # Mise à jour des collisions avec les murs de la boîte
        for i in range(Natoms):
            loc = apos[i]
            if abs(loc.x) > L / 2:
                if loc.x < 0:
                    p[i].x = abs(p[i].x)
                else:
                    p[i].x = -abs(p[i].x)
            if abs(loc.y) > L / 2:
                if loc.y < 0:
                    p[i].y = abs(p[i].y)
                else:
                    p[i].y = -abs(p[i].y)
        print(4)
        # Mise à jour des collisions entre les sphères
        for ij in hitlist:
            i, j = ij
            ptot = p[i] + p[j]
            mtot = 2 * mass
            Vcom = ptot / mtot
            posi = apos[i]
            posj = apos[j]
            vi = p[i] / mass
            vj = p[j] / mass
            rrel = posi - posj
            vrel = vj - vi
            print(5)

            if vrel.mag2 == 0 or rrel.mag > Ratom:
                continue

            dx = dot(rrel, vrel.hat)
            dy = cross(rrel, vrel.hat).mag
            alpha = asin(dy / (2 * Ratom))
            d = (2 * Ratom) * cos(alpha) - dx
            deltat = d / vrel.mag

            posi = posi - vi * deltat
            posj = posj - vj * deltat
            pcomi = p[i] - mass * Vcom
            pcomj = p[j] - mass * Vcom
            rrel = hat(rrel)
            pcomi = pcomi - 2 * dot(pcomi, rrel) * rrel
            pcomj = pcomj - 2 * dot(pcomj, rrel) * rrel
            p[i] = pcomi + mass * Vcom
            p[j] = pcomj + mass * Vcom
            apos[i] = posi + (p[i] / mass) * deltat
            apos[j] = posj + (p[j] / mass) * deltat
            print(6)
        simulation_time += dt

    # Convertir les listes en tableaux NumPy
    collision_times_array = np.array(collision_times)
    collision_distances_array = np.array(collision_distances)

    return collision_times_array, collision_distances_array

# Définir la variable running avant d'appeler la fonction
running = True

# Appel de la fonction et impression des données
collision_times_array, collision_distances_array = follow_the_bouncing_ball()
print("Collision Times:", collision_times_array)
print("Collision Distances:", collision_distances_array)

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre


def stop_simulation(evt):
    global running
    if evt.key == "esc":  # Vérifie si la touche est 'Esc'
        running = False  # Arrête la simulation


# Désactive la saisie de texte dans la fenêtre d'animation
#def disable_default_text_input(evt):
#    pass  # Ne fait rien, empêche l'entrée de texte par défaut


# Liaison des touches pour empêcher les inputs
#animation.bind(
#    "keydown", disable_default_text_input
#)  # Empêche l'affichage de caractères
#animation.bind("keydown", stop_simulation)  # Détecte `Esc`

running = True
print("Collision Times 2:")
print("Collision Distances 2:")
while running:
    rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

    frame_counter = 0

    if frame_counter % 1000 == 0:  # Enregistre toutes les 1000 itérations
        save_quantite_mouvement()

    frame_counter += 1



    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    vitesse = []  # vitesse instantanée de chaque sphère
    deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for i in range(Natoms):
        vitesse.append(
            p[i] / mass
        )  # par définition de la quantité de nouvement pour chaque sphère
        deltax.append(
            vitesse[i] * dt
        )  # différence avant pour calculer l'incrément de position
        Atoms[i].pos = apos[i] = (
            apos[i] + deltax[i]
        )  # nouvelle position de l'atome après l'incrément de temps dt

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE ####
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L / 2:
            if loc.x < 0:
                p[i].x = abs(p[i].x)  # renverse composante x au mur de gauche
            else:
                p[i].x = -abs(p[i].x)  # renverse composante x au mur de droite
        if abs(loc.y) > L / 2:
            if loc.y < 0:
                p[i].y = abs(p[i].y)  # renverse composante y au mur du bas
            else:
                p[i].y = -abs(p[i].y)  # renverse composante y au mur du haut

    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist = checkCollisions()

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for ij in hitlist:
        # définition de nouvelles variables pour chaque paire de sphères en collision
        i = ij[0]  # extraction du numéro des 2 sphères impliquées à cette itération
        j = ij[1]
        ptot = p[i] + p[j]  # quantité de mouvement totale des 2 sphères
        mtot = 2 * mass  # masse totale des 2 sphères
        Vcom = (
            ptot / mtot
        )  # vitesse du référentiel barycentrique/center-of-momentum (com) frame
        posi = apos[i]  # position de chacune des 2 sphères
        posj = apos[j]
        vi = p[i] / mass  # vitesse de chacune des 2 sphères
        vj = p[j] / mass
        rrel = posi - posj  # vecteur pour la distance entre les centres des 2 sphères
        vrel = vj - vi  # vecteur pour la différence de vitesse entre les 2 sphères

        # exclusion de cas où il n'y a pas de changements à faire
        if vrel.mag2 == 0:
            continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        if rrel.mag > Ratom:
            continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        dx = dot(
            rrel, vrel.hat
        )  # rrel.mag*cos(theta) où theta is the angle between vrel and rrel:
        dy = cross(rrel, vrel.hat).mag  # rrel.mag*sin(theta)
        alpha = asin(
            dy / (2 * Ratom)
        )  # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        d = (2 * Ratom) * cos(
            alpha
        ) - dx  # distance traveled into the atom from first contact
        deltat = (
            d / vrel.mag
        )  # time spent moving from first contact to position inside atom

        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        posi = posi - vi * deltat  # back up to contact configuration
        posj = posj - vj * deltat
        pcomi = (
            p[i] - mass * Vcom
        )  # transform momenta to center-of-momentum (com) frame
        pcomj = p[j] - mass * Vcom
        rrel = hat(rrel)  # vecteur unitaire aligné avec rrel
        pcomi = (
            pcomi - 2 * dot(pcomi, rrel) * rrel
        )  # bounce in center-of-momentum (com) frame
        pcomj = pcomj - 2 * dot(pcomj, rrel) * rrel
        p[i] = pcomi + mass * Vcom  # transform momenta back to lab frame
        p[j] = pcomj + mass * Vcom
        apos[i] = (
            posi + (p[i] / mass) * deltat
        )  # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération
        apos[j] = posj + (p[j] / mass) * deltat
