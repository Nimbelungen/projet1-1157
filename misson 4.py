import math

masse_barge = 10
masse_grue = 15
masse_poid_deplace = 1.5
distance_deplacement = 2
hauteur_deplacement = 0
g = 9.81
longeur_barge = 1
largeur_barge = 1
hauteur_barge = 0.8
hauteur_grue = 0.5

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
    """
    Calculer la hauteur immergée de la barge
    :return: La hauteure hc immergée de la barge
    """
    v_immerge = ((masse_barge + masse_grue) * g) / (1000 * g)
    hc = v_immerge / (longeur_barge * largeur_barge)
    return hc


def angle_max():
    """
    :return: l'inclinaison maximum avant que la barge ne chavire à coup sur
    """
    thetamax = math.atan((hauteur_barge - calcul_hc()) / (longeur_barge / 2))
    return thetamax


def center_gravity():
    """
    :return: Retoune les coordonnées du centre de gravité
    """
    center_gravity_x = ((masse_barge * 0) + (masse_grue * 0) + (masse_poid_deplace * distance_deplacement)) / (masse_barge + masse_grue + masse_poid_deplace)
    center_gravity_y = ((masse_barge * 0) + (masse_grue * 0) + (masse_poid_deplace * 0)) / (masse_barge + masse_grue + masse_poid_deplace)
    center_gravity_z = ((masse_barge * (- hauteur_barge / 2)) + (masse_grue * (hauteur_barge - calcul_hc() + hauteur_grue)) + (masse_poid_deplace * hauteur_deplacement)) / (masse_barge + masse_grue + masse_poid_deplace)
    return [center_gravity_x, center_gravity_z, center_gravity_y]


def calcul_theta():
    # 1 : Calcul du nouveau centre de gravité
    # 2 : Choix de theta
    # 3 : Calcul du nouveau
    for i in range(1000000):
        ca = masse_poid_deplace * g * distance_deplacement
        cr = math.sqrt()
        if ca == cr:
            return


print(calcul_hc())
