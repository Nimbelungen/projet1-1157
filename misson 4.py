import math

# Variables a definir
masse_barge = 10
masse_grue = 15
masse_poid_deplace = 1.5
distance_deplacement = 2
hauteur_deplacement = 0
longeur_barge = 1
largeur_barge = 1
hauteur_barge = 0.8
hauteur_grue = 0.5

# Variables fixes
g = 9.81
tan_theta = 0

""" Initialiser les variables avec les input
masse_barge = float(input("rentrer le masse de la barge :"))
masse_grue = float(input("rentrer la masse de la grue :"))
masseC = float(input("rentrer la masse du poid déplacé:"))
d = float(input("la distance entre l'objet déplacé et la barge :"))
g = 9.81
l1 = float(input("longueur de la base :" ))
l2 = float(input("largeur de la base : "))
h1 = float(input("hauteur de toute la barge"))
"""

"""ATTENTION : L ENSEMBLE DES CALCULS CE FAIT DANS LE REPERE OU
L AXE X = parrallèle à la surface de l'eau, 0 au milieu de la barge (et de la grue)
L AXE Y = perpendiculaire à x (sur l'eau), 0 au milieu de la barge (et de la grue)
L AXE Z = verticale (perendiculaire à l'eau), 0 au niveau de l'eau
"""


def calcul_hc():
    """Calculer la longeur immergee
    :return: la distence hc, ou hc est la longeur du volume immerge
    """
    sub_volume = (masse_grue + masse_barge + masse_poid_deplace) / 1000
    hc = sub_volume / (longeur_barge * largeur_barge)
    return hc


def angle_max(axe):
    """
    Calcul l'angle d'inclinaison maximum de la barge. Si cette valeur est dépassée, l'au s'inflitre sur la barge
    :type axe: str
    :param axe: l'axe autaur duquel la barge va trouner : si toune autour de y => longeur ; si tourne autour de x => largeur
    :return: l'angle maximum EN RADIANS
    """
    global tan_theta
    if axe == "y":
        tan_theta = (hauteur_barge - calcul_hc()) / (longeur_barge / 2)
    elif axe == "x":
        tan_theta = (hauteur_barge - calcul_hc()) / (largeur_barge / 2)
    # elif axe == "xy": todo trouver comment faire si ça tourne autour de 2 axes => si la grua à une base qui bouge

    theta = math.atan2(tan_theta)
    return theta


def center_gravity_global(m1, m2, m3, d1, d2, d3):
    """
    :param m1:
    :param m2:
    :param m3:
    :param d1: tuples (vecteur depuis l'origine) (coordonné du centre de gravité de m1)
    :param d2: tuples (vecteur depuis l'origine) (coordonné du centre de gravité de m2)
    :param d3: tuples (vecteur depuis l'origine) (coordonné du centre de gravité de m3)
    :return:
    """
    cg = []
    for axes in range(2):
        coord = (m1 * d1[axes]) + (m2 * d2[axes]) + (m3 * d3[axes]) / m1 + m2 + m3
        cg.append(coord)
    return tuple(cg)


def center_gravity(initial):
    """
    Calcule les coordonnées du centre de gravité (voir repère)
    :type initial: bool
    :param initial: True => situation initiale / False => situation avec dé^placement
    :return: Les cooronées du centre de gravite en fonction de initial
    """



def center_pousse_tra(parallel1, parallel2, hauteur, theta):
    code



def calcul_theta():
    # 1 : Calcul du nouveau centre de gravite
    # 2 : Choix de theta
    # 3 : Calcul du nouveau centre de gravite
    # 4 : Calcul du nouveau centre de pousse
    code
