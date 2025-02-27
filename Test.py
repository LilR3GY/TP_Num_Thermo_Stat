from vpython import *
import numpy as np
from scipy.stats import maxwell
import matplotlib.pyplot as plt  # Pour le tracé du graphique

# Déclaration de variables influençant le temps d'exécution de la simulation
Nions = 25  # Nombre de cœurs ioniques
Rion = 0.02  # Taille fictive des ions
Natoms = 200  # Nombre d'atomes
dt = 1E-8  # Pas d'incrémentation temporelle

# Déclaration de variables physiques
DIM = 2  # Nombre de degrés de liberté de la simulation
mass = 9.109e-31  # Masse d'un atome
Ratom = 0.01  # Taille des atomes
k = 1.4E-23  # Constante de Boltzmann
T = 1000  # Température ambiante augmentée à 1000K

#### CANEVAS DE FOND ####
L = 1  # Conteneur en forme de cube
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
    Atoms.append(simple_sphere(pos=vector(x, y, z), radius=Ratom, color=gray))
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

# Fonction pour identifier les collisions entre électrons et ions
def checkIonCollisions():
    hitlist = []
    r2 = (Ratom + Rion) ** 2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(len(Ions)):
            aj = ion_pos[j]
            dr = ai - aj
            if mag2(dr) < r2:
                hitlist.append([i, j])
    return hitlist

# Fonction pour générer une quantité de mouvement aléatoire
def generateRandomMomentum():
    speed = maxwell.rvs(scale=np.sqrt(k * T / mass))
    theta = 2 * np.pi * random()
    px = mass * speed * np.cos(theta)
    py = mass * speed * np.sin(theta)
    return vector(px, py, 0)

# Listes pour stocker les données
avg_momentum_list = []
single_electron_momentum = []
avg_direction_list = []  # Stocke la direction moyenne des électrons
single_electron_direction = []  # Stocke la direction d'un seul électron

# Boucle principale
Valeurs = 500

for iteration in range(Valeurs):
    rate(1000)

    # Déplace les sphères
    vitesse = []
    deltax = []
    for i in range(Natoms):
        vitesse.append(p[i] / mass)
        deltax.append(vitesse[i] * dt)
        Atoms[i].pos = apos[i] = apos[i] + deltax[i]

    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L / 2:
            p[i].x *= -1
            apos[i].x = L / 2 * (1 if loc.x > 0 else -1)
        if abs(loc.y) > L / 2:
            p[i].y *= -1
            apos[i].y = L / 2 * (1 if loc.y > 0 else -1)

    # Identification des collisions avec les ions
    hitlist = checkIonCollisions()

    # Gestion des collisions
    for ij in hitlist:
        i = ij[0]
        j = ij[1]
        p[i] = generateRandomMomentum()
        apos[i] += (p[i] / mass) * dt

    # Calcul de la quantité de mouvement moyenne
    avg_momentum = np.mean([mag(p_i) for p_i in p])
    avg_momentum_list.append(avg_momentum)

    # Stocke la quantité de mouvement du premier électron
    single_electron_momentum.append(mag(p[0]))

    # Calcul de la direction moyenne des électrons
    avg_direction = np.mean([np.arctan2(p_i.y, p_i.x) for p_i in p])
    avg_direction_list.append(avg_direction)

    # Stocke la direction de l'électron choisi
    single_electron_direction.append(np.arctan2(p[0].y, p[0].x))

    # Recalcul de la température toutes les 10 itérations
    if iteration % 10 == 0:
        if avg_momentum > 1E-32:  # Si la quantité de mouvement moyenne est significativement non nulle
            # Calcul de la nouvelle température basée sur la quantité de mouvement
            T = (avg_momentum ** 2) * mass / (k * (DIM / 2))
        else:
            print("La quantité de mouvement est presque nulle, arrêt de la simulation.")
            break


def plot_momentum_decay(p0, tau, t_max, num_points=1000):
    t = np.linspace(0, t_max, num_points)  # Génère des points de temps
    p_t = p0 * np.exp(-t / tau)  # Calcule <p(t)>

    # Tracer la fonction de décroissance
    plt.plot(t, p_t, label='Décroissance Exponentielle', color='blue')

# Paramètres pour le tracé
p0 = avg_momentum_list[0] if avg_momentum_list else 1.0  # Valeur initiale de <p(t0)>
tau = 0.00000000001
t_max = Valeurs / 5000000
num_points = 100  # Nombre de points pour le tracé

start = 0  # Début
stop = Valeurs / 10  # Fin
step = 0.1  # Pas

# Calculer le nombre de points nécessaires
num_points = int((stop - start) / step)

# Utiliser linspace pour générer la liste
values = np.linspace(start, stop, num_points)

plt.figure() 
# Appel de la fonction pour tracer la décroissance
plot_momentum_decay(p0, tau, t_max, num_points)

# Tracé du graphique de la quantité de mouvement moyenne en fonction du temps
time_values = [dt * i for i in range(len(avg_momentum_list))]
plt.plot(time_values, avg_momentum_list)
plt.xlabel('Temps (s)')
plt.ylabel('Quantité de mouvement moyenne (kg.m/s)')
plt.title('Quantité de mouvement moyenne des électrons en fonction du temps')

plt.grid()
plt.show()