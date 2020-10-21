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
lb = 0      # longeur_barge
hb = 0      # hauteur_barge
lg = 0      # longeur_grue
hg = 0      # hauteur_grue
dd = 0      # distance_deplacement
theta = 0   # Angle d'inclinaison de la barge
# --- Variables immuables ---
g = 9.81    # Accélération gravitationelle

