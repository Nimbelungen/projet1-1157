import math

"""ATTENTION : L ENSEMBLE DES CALCULS CE FAIT DANS LE REPERE OU
L AXE X = parrallèle à la surface de l'eau, 0 au milieu de la barge (et de la grue)
L AXE Z = verticale (perendiculaire à l'eau), 0 au niveau de l'eau
"""
print("On considère que la masse de la grue a pour position x d/2 après déplacement du poid")

# ------ Definir les variables ------
# --- Variables muablbes ---
mb = 0      # masse_barge
mg = 0      # masse_grue
mp = 0      # masse_poid
# mcp = 0     # masse_contre_poid
lb = 0      # longeur_barge - lorsqu un volume est necessaire, on considere la longeur = largeur
hb = 0      # hauteur_barge
lg = 0      # longeur_grue
hg = 0      # hauteur_grue
dd = 0  # distance_deplacement
hp = 0  # Hauteur du poid lors de la situation initiale (il sera en z = hd après le déplacement)
theta = 0  # Angle d'inclinaison de la barge
# --- Variables immuables ---
g = 9.81    # Accélération gravitationelle


# ------ Les fonctions elementaires ------
def find_hc():
    """
    Calculate the submerged length of the barge
    :return: If hc < barge height : the distance hc, where hc is the length follow the submerged z-axis of the barge. Otherwise, False
    """
    sub_volume = (mg + mb + mp) / 1000
    hc = sub_volume / (2 * lb)
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
    angle_max1 = math.atan2(tan_theta1)
    tan_theta2 = find_hc() / (lb / 2)
    angle_max2 = math.atan2(tan_theta2)
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
    for axes in range(1):
        coord = (m1 * d1[axes]) + (m2 * d2[axes]) + (m3 * d3[axes]) / m1 + m2 + m3
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
        d3 = (0, hd + hp)
    else:
        d2 = (dd / 2, hd + hg)
        d3 = (dd, hd)
    return center_gravity_glob(mb, mg, mp, d1, d2, d3)


def center_thrust(init, parrallel_left, parrallel_right, height, angle):
    """ Caution : Rotation is also applied to the axes!) - CENTRE DE POUSSEE
    :type init: bool
    :type parrallel_left: float
    :type parrallel_right: float
    :type height: float
    :type angle: float
    :return:
    """
    hc = find_hc()
    if init:
        ctx = lb / 2
        ctz = hc / 2
        return tuple([ctx, ctz])
    else:
        # For more information, see the README file
        dist_pr = (height / 3) * ((parrallel_right + (2 * parrallel_left)) / (
                    parrallel_right + parrallel_left))  # Formula find on Wikipedia
        c1 = [-lb / 2, -hc]
        c2 = [lb / 2, -hc]
        d1 = (hc - (math.tan(angle) * lb / 2)) / 2
        d2 = (hc + (math.tan(angle) * lb / 2)) / 2
        p1 = [c1[0], c1[1] + d1]
        p2 = [c2[0], c2[1] + d2]
        # droite_p1_p2 = z - p1[1] = ( (p2[1] - p1[1]) / (p2[0] - p1[0]) ) * (x - p1[0])
        # droite_perp_dist_pr = z = dist_pr
        # On trouve donc les coordonnées qui sont [x, z] telle que ces 2 équations soient vérifiées
        x_center_thust = (lb / 2) - dist_pr
        z_center_thust = (((p2[0] - p1[0]) / (p2[1] - p1[1])) * (dist_pr - p1[1])) - p1[0]
        return tuple([x_center_thust, z_center_thust])


def center_thrust_theta(angle):
    """
    :type angle: float
    :return:
    """
    init = False
    # find parrallel_left
    # find parrallel_right
    # find height

    # Coordonnes dans le repere non 'vertical'
    coordonate_t = center_thrust(False, x, x, x, angle)
    # faire tourner le repere
    coordonate_l = list(coordonate_t)

    return coordonate_l
