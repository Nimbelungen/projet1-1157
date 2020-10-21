
import math

masse1 = float(input ("rentrer le masse de la barge :"))
masse2 = float(input ("rentrer la masse de la grue :"))
masseC = float(input ("rentrer la masse du poid déplacé:"))
d = float(input ("la distance entre l'objet déplacé et la barge :"))
g = 9.81
l1 = float(input ("longueur de la base :" ))
l2 = float(input ("largeur de la base : "))
h1 = float(input ("hauteur de toute la barge"))

def calcul_hc():
    V = ((masse1+masse2)*g)/(1000*g)
    hc =V/(l1*l2)
    return hc

def angle_max():
    thetamax = math.atan((h1-calcul_hc())/(l1/2))
    return thetamax


def calcul_theta():
    for i in range(1000000):
        ca = masseC*g*d
        cr = math.sqrt()
        if ca == cr :
            return 


print(calcul_hc())
