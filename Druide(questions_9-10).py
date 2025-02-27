from vpython import *
import numpy as np

# Déclaration de variables influençant le temps d'exécution de la simulation
Nions = 25  # Nombre de cœurs ioniques
Rion = 0.02  # Taille fictive des ions
Natoms = 200  # Change this to have more or fewer atoms
dt = 1E-5  # Pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
DIM = 2  # Nombre de degrés de liberté de la simulation 
mass = 4E-3 / 6E23  # Masse d'un atome (remplacé par un autre type si nécessaire)
Ratom = 0.01  # Taille des atomes
k = 1.4E-23  # Constante de Boltzmann
T = 300  # Température ambiante
r2 = (Ratom + Rion) ** 2  # Distance critique pour collision

# Déclaration des variables du champ électrique
E_intensity = -1E-29  # Intensité du champ électrique (modifiable dynamiquement)
E_direction = 'vertical'  # 'horizontal' ou 'vertical'
charge_electron = -1.6E-19  # Charge d'un électron (Coulombs)

#### CANEVAS DE FOND ####
L = 1  # Conteneur en forme de cube L sur un côté
gray = color.gray(0.7)  # Couleur des bords du conteneur
animation = canvas(width=750, height=500)
animation.range = L

#### ARÊTES DE BOÎTE 2D ####
d = L / 2 + Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d, -d, 0), vector(d, -d, 0), vector(d, d, 0), vector(-d, d, 0), vector(-d, -d, 0)])

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Atoms = []
p = []
apos = []
pavg = sqrt(2 * mass * (DIM / 2) * k * T)

for i in range(Natoms):
    x = L * random() - L / 2
    y = L * random() - L / 2
    z = 0
    Atoms.append(sphere(pos=vector(x, y, z), radius=Ratom, color=gray)) 
    apos.append(vec(x, y, z)) 
    phi = 2 * pi * random()  
    px = pavg * cos(phi)  
    py = pavg * sin(phi)
    p.append(vector(px, py, 0)) 

# Initialisation des cœurs ioniques
Ions = []
ion_pos = []

grid_size = int(sqrt(Nions))
dist = L / grid_size

for i in range(grid_size):
    for j in range(grid_size):
        if len(Ions) >= Nions:
            break
        x = -L / 2 + i * dist + dist / 2
        y = -L / 2 + j * dist + dist / 2
        ion = simple_sphere(pos=vector(x, y, 0), radius=Rion, color=color.red)
        Ions.append(ion)
        ion_pos.append(vector(x, y, 0))

# Fonction pour identifier les collisions entre les électrons et les cœurs ioniques
def checkIonCollisions():
    hitlist = []  # Initialisation
    for i in range(Natoms):
        ai = apos[i]
        for j in range(len(Ions)):
            aj = ion_pos[j]
            dr = ai - aj  # Distance vectorielle
            if mag2(dr) < r2:  # Test de collision
                hitlist.append([i, j])  # Liste des paires en collision
    return hitlist

# Fonction pour générer une quantité de mouvement aléatoire selon Maxwell-Boltzmann
def generateRandomMomentum():
    # Direction aléatoire
    theta = np.arccos(2 * random() - 1)  # Azimuthal angle
    phi = 2 * np.pi * random()  # Polar angle
    # Norme de la quantité de mouvement basée sur la distribution de Maxwell-Boltzmann
    momentum_magnitude = np.sqrt(-2 * mass * k * T * np.log(1 - random()))  # Distribution exponentielle
    px = momentum_magnitude * np.sin(theta) * np.cos(phi)
    py = momentum_magnitude * np.sin(theta) * np.sin(phi)
    return vector(px, py, 0)

# Fonction pour appliquer le champ électrique
def apply_electric_field(p, dt, hitlist, Natoms, mass):
    global E_intensity  # Permettre la modification dynamique de l'intensité
    E_field = vector(E_intensity, 0, 0) if E_direction == 'horizontal' else vector(0, E_intensity, 0)
    hit_indices = {collision[0] for collision in hitlist}  # Utilisation d'un ensemble pour optimisation
    for i in range(Natoms):
        if i not in hit_indices:  # Vérifier les atomes hors collision
            F_electric = charge_electron * E_field  # Force de Coulomb
            p[i] += F_electric * dt / mass  # Modification de la quantité de mouvement

# Variables pour stocker les positions moyennes au fil du temps
positions_paralleles = []
positions_perpendiculaires = []

# Fonction pour calculer la position moyenne des électrons
def compute_average_position():
    avg_pos_parallel = 0
    avg_pos_perpendicular = 0
    for i in range(Natoms):
        # Si champ horizontal, la composante parallèle est la position x
        if E_direction == 'horizontal':
            avg_pos_parallel += apos[i].x
            avg_pos_perpendicular += apos[i].y
        else:  # Si champ vertical, la composante parallèle est la position y
            avg_pos_parallel += apos[i].y
            avg_pos_perpendicular += apos[i].x
    
    avg_pos_parallel /= Natoms
    avg_pos_perpendicular /= Natoms
    
    return avg_pos_parallel, avg_pos_perpendicular

# Ajout dans la boucle principale
dt = 1E-5  # Pas d'incrémentation temporelle
max_iterations = 1000  # Limite d'itérations pour éviter une boucle infinie
iteration = 0

while iteration < max_iterations:
    rate(500)  # Réduction de la vitesse de simulation pour alléger la charge CPU

    # Déplacement des sphères
    for i in range(Natoms):
        vitesse = p[i] / mass
        deltax = vitesse * dt
        Atoms[i].pos = apos[i] = apos[i] + deltax

    # Gestion des collisions avec les parois
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L / 2:
            p[i].x *= -1
            apos[i].x = (L / 2) * (1 if loc.x > 0 else -1)
        if abs(loc.y) > L / 2:
            p[i].y *= -1
            apos[i].y = (L / 2) * (1 if loc.y > 0 else -1)

    # Identification des collisions avec les cœurs ioniques
    hitlist = checkIonCollisions()

    # Application du champ électrique sur les atomes sans collision
    apply_electric_field(p, dt, hitlist, Natoms, mass)
    
    # Gestion des collisions inélastiques
    for ij in hitlist:
        i = ij[0]
        p[i] = 0.5 * p[i] + 0.5 * generateRandomMomentum()

    # Calcul des positions moyennes
    avg_pos_parallel, avg_pos_perpendicular = compute_average_position()
    positions_paralleles.append(avg_pos_parallel)
    positions_perpendiculaires.append(avg_pos_perpendicular)
    
    iteration += 1

    import matplotlib.pyplot as plt

# Générer l'axe du temps
time = np.arange(max_iterations) * dt

# Tracer les positions moyennes parallèle et perpendiculaire au champ
plt.figure(figsize=(10, 5))

# Composante parallèle
plt.subplot(1, 2, 1)
plt.plot(time, positions_paralleles, label='Composante parallèle')
plt.xlabel('Temps (s)')
plt.ylabel('Position moyenne')
plt.title('Composante parallèle au champ électrique')
plt.legend()

# Ajouter une flèche pour le sens du champ
if E_direction == 'horizontal':
    if E_intensity > 0:
        arrow(pos=vector(0, 0.6, 0), axis=vector(0.5 * L, 0, 0), shaftwidth=0.05 * L, color=color.red)
    else:
        arrow(pos=vector(0, 0.6, 0), axis=vector(-0.5 * L, 0, 0), shaftwidth=0.05 * L, color=color.red)
    label = text(text='Champ horizontal', pos=vector(-0.25 * L, 0.68 * L, 0), height=0.05, color=color.red)
else:  # Cas pour un champ vertical
    if E_intensity > 0:
        arrow(pos=vector(-0.6, 0, 0), axis=vector(0, 0.5 * L, 0), shaftwidth=0.05 * L, color=color.red)
    else:
        arrow(pos=vector(-0.6, 0, 0), axis=vector(0, -0.5 * L, 0), shaftwidth=0.05 * L, color=color.red)
    label = text(text='Champ vertical', pos=vector(-0.68 * L, -0.20 * L, 0), height=0.05, color=color.red, axis=vector(0, 1, 0))

# Composante perpendiculaire
plt.subplot(1, 2, 2)
plt.plot(time, positions_perpendiculaires, label='Composante perpendiculaire')
plt.xlabel('Temps (s)')
plt.ylabel('Position moyenne')
plt.title('Composante perpendiculaire au champ électrique')
plt.legend()

# Affichage du graphique
plt.tight_layout()
plt.show()