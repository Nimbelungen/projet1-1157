import math

"""ATTENTION : L ENSEMBLE DES CALCULS CE FAIT DANS LE REPERE OU
L AXE X = parrallèle à la surface de l'eau, 0 au milieu de la barge (et de la grue)
L AXE Z = verticale (perendiculaire à l'eau), 0 au niveau de l'eau
"""

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
dd = 0      # distance_deplacement
theta = 0   # Angle d'inclinaison de la barge
# --- Variables immuables ---
g = 9.81    # Accélération gravitationelle


# ------ Les fonctions elementaires ------
def return_hc():
    """Calculer la longeur immergee de la barge
        :return: Si hc < hauteur barge : la distence hc, ou hc est la longeur suivent l'axe z immergee de la barge
                    Sinon, False
        """
    sub_volume = (mg + mb + mp) / 1000
    hc = sub_volume / (2 * lb)
    if hc < hb:
        return hc
    else:
        return False


def return_angle_max():
    """
    Calcul l'angle d'inclinaison maximum de la barge. Si cette valeur est dépassée, l'au s'inflitre sur la barge
    :return: l'angle maximum EN RADIANS (selon l axe x - cad rotation autour de l axe y, imaginaire en 2D)
    """
    tan_theta = (hb - return_hc()) / (lb / 2)
    angle_max = math.atan2(tan_theta)
    return angle_max


