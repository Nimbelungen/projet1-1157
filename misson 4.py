import math

masse_barge = 10
masse_grue = 15
masse_poid_deplace = 1.5
distance_deplacement = 2
g = 9.81
longeur_barge = 1
largeur_barge = 1
hauteur_barge = 0.8

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


def calcul_hc():
    V = ((masse_barge + masse_grue) * g) / (1000 * g)
    hc =V/(longeur_barge * largeur_barge)
    return hc


def angle_max():
    thetamax = math.atan((hauteur_barge - calcul_hc()) / (longeur_barge / 2))
    return thetamax


def calcul_theta():
    for i in range(1000000):
        ca = masse_poid_deplace * g * distance_deplacement
        cr = math.sqrt()
        if ca == cr:
            return


print(calcul_hc())
