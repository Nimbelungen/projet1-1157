import math
""" INFORMATIONS
Ceci est l'exercice de phyqique H2. La méthode de résolution de l'exercice est la suivente:
1. La masse mp n'est pas considérée comme négligeable
2. Le déplacement de mp provoque un changement de position du centre de gravité
3. On considère alors que la barge va subir une inclinaison d'angle thêta
4. Le centre de poussé change alors aussi d'emplacement
5. La barge sera à l'équilibre lorsque les positions en x du centre de gravité et du centre de poussée seront =
6. On trouve l'angle par essaie-erreur car les coordonnées du centre de poussées dépendent direcetement de cet angle. 
"""

"""ATTENTION : L ENSEMBLE DES CALCULS SE FONT DANS LE REPERE OU:
L AXE X = parrallèle à la surface de l'eau, 0 au milieu de la barge (et de la grue)
L AXE Z = verticale (perendiculaire à l'eau), 0 au niveau de l'eau
(sauf pour la fonction 'center_thust()', voir axes verts dans le README)
"""
print("On considere que la masse de la grue a pour position x d/2 après deplacement du poids")

# ------ Definir les variables ------
# --- Variables muablbes ---
mb = 20.0  # masse_barge
mg = 13.5  # masse_grue
mp = 1.5  # masse_poid
# mcp = 0   # masse_contre_poid
lb = 1.0  # longeur_barge - lorsqu un volume est necessaire, on considere la longeur = largeur
hb = 0.08  # hauteur_barge
lg = 0.0  # longeur_grue
hg = 0.1  # hauteur_grue
dd = 2.0  # distance_deplacement
hp = 0  # Hauteur du poid lors de la situation initiale (il sera en z = hd après le déplacement)
theta = 0  # Angle d'inclinaison de la barge
# --- Variables immuables ---
g = 9.81  # Accélération gravitationelle


# ------ Les fonctions elementaires ------
def find_hc():
    """
    Calculate the submerged length of the barge
    :return: If hc < barge height : the distance hc, where hc is the length follow the submerged z-axis of the barge. Otherwise, False
    """
    sub_volume = (mg + mb + mp) / 1000
    hc = sub_volume / (lb ** 2)
    if hc < hb:
        return hc
    else:
        return False


def angle_max():
    """
    Calculate the maximum tilt angle of the barge. If this value is exceeded, water will seep onto the barge.
    :return: the maximum angle IN RADIANS (along x-axis - cad rotation around y-axis, imaginary in 2D)
    """
    # There is 2 critical points, the critical point is the minimum
    tan_theta1 = (hb - find_hc()) / (lb / 2)
    angle_max1 = math.atan(tan_theta1)
    tan_theta2 = find_hc() / (lb / 2)
    angle_max2 = math.atan(tan_theta2)
    if angle_max1 <= angle_max2:
        return angle_max1
    else:
        return angle_max2


def center_gravity_glob(m1, m2, m3, d1, d2, d3):
    """
    :type m1: float
    :type m2: float
    :type m3: float
    :type d1: tuple
    :type d2: tuple
    :type d3: tuple
    :return: a tuple with the coordinates of the center of gravity og m1, m2 and m3
    """
    cg = []
    for axes in [0, 1]:
        coord = ((m1 * d1[axes]) + (m2 * d2[axes]) + (m3 * d3[axes])) / (m1 + m2 + m3)
        cg.append(coord)
    return tuple(cg)


def center_gravity(init):
    hc = find_hc()
    """
    :type init: bool
    :return: tuple with the coordinates of the center of gravity og mb, mg and mp depending on the situation
    """
    dist_axex_cgbarge = (hb / 2) - hc
    hd = hb - hc
    d1 = (0, dist_axex_cgbarge)
    if init:
        d2 = (0, hd + hg)
        d3 = (0, hd + hg)
    else:
        d2 = (0, hd + hg)
        d3 = (dd, hd + hg)
    return center_gravity_glob(mb, mg, mp, d1, d2, d3)


def center_thrust(init, angle):
    """ Caution : Rotation is also applied to the axes!) - CENTRE DE POUSSEE
    :type init: bool
    :type angle: float
    :return: Un tuple avec les coordonnees (x, z) du centre de poussee dans un repere tourne vers la droite d un angle
    'angle'
    """
    hc = find_hc()
    if init:
        ctx = lb / 2
        ctz = hc / 2
        return tuple([ctx, ctz])
    else:
        # variables utiles - For more information, see the README file
        parrallel_left = (hc - (math.tan(angle) * (lb / 2)))  # longeur
        parrallel_right = (hc + (math.tan(angle) * (lb / 2)))  # longeur

        c1 = [-lb / 2, -hc]
        c2 = [lb / 2, -hc]
        d1 = parrallel_left / 2
        d2 = parrallel_right / 2
        p1 = [c1[0], c1[1] + d1]
        p2 = [c2[0], c2[1] + d2]

        height = lb  # math.sqrt(((c2[0] - c1[0]) ** 2) + ((c2[1] - c2[1]) ** 2))

        dist_pr = (height / 3) * ((parrallel_right + (2 * parrallel_left)) / (parrallel_right + parrallel_left))  #
        # Formula find on Wikipedia

        # droite_p1_p2 = z - p1[1] = ( (p2[1] - p1[1]) / (p2[0] - p1[0]) ) * (x - p1[0])
        # droite_perp_dist_pr = z = dist_pr
        # On trouve donc les coordonnées qui sont [x, z] telle que ces 2 équations soient vérifiées
        x_center_thust = (lb / 2) - dist_pr
        z_center_thust = (((p2[1] - p1[1]) / (p2[0] - p1[0])) * (x_center_thust - p1[0])) + p1[1]
        return tuple([x_center_thust, z_center_thust])


def rotate_center_thust(angle):
    """
    :type angle: float
    :return: un tuple (x, y) avec les coordonnee du centre de poussee dans le repere initiale
    """
    init = False
    # Coordonnes dans le repere non 'vertical'
    coordonate_t = center_thrust(False, angle)
    # faire tourner le repere
    coordonate_l = tuple([(coordonate_t[0] * math.cos(-angle)) - (coordonate_t[1] * math.sin(-angle)),
                          (coordonate_t[0] * math.sin(-angle)) + (coordonate_t[1] * math.cos(-angle))])
    return coordonate_l


# ------ Trouver l'angle ------
def find_theta():
    """
    Algorithme de recherche de l'angle theta
    :return: l'angle theta en radian
    """
    first = 0.0
    last = angle_max()
    find = False

    while first <= last and not find:
        middle = (first + last) / 2
        if rotate_center_thust(middle)[0] == center_gravity(False)[0]:
            return middle
        else:
            if rotate_center_thust(middle)[0] < center_gravity(False)[0]:
                first = middle + 0.00000000000000001

            else:
                last = middle - 0.00000000000000001


def to_degrees(angle):
    return (180 * angle) / math.pi


# ------ Test ------
print("Hc :", find_hc())
print("Angle max :", angle_max())
print("Centre de gravité initiale :", center_gravity(True))
print("Centre de gravité après mouvement :", center_gravity(False))
print("Centre de poussé après mouvement : ", rotate_center_thust(0.03697))
print("L'angle theta vaut :", find_theta(), "rad")
print("L'angle theta vaut :", to_degrees(find_theta()), "°")

print("")
print("----- Corrections -----")
print("Cg x", center_gravity(False)[0])
print("Cp x", rotate_center_thust(find_theta())[0])
