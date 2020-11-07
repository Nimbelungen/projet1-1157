from formulas import *

# ---- Variables ----
# <name>_x = length [m]
# <name>_y = width  [m]
# <name>_z = height [m]
# <name>_mass = mass    [kg]

# -- Barge --
barge_x = 0.5  # Barge length
barge_y = 0.5  # Barge width
barge_z = 0.3  # Barge height
barge_mass = 0.3  # Barge mass

# -- Grue --
# First Piece
grue1_x = 0.15  # Grue First Piece length
grue1_z = 0.05  # Grue First Piece height
grue1_mass = 0.05  # Grue First Piece mass
# Grue First Piece angle is always 0
# Second Piece
grue2_x = 0.350  # Grue Second Piece length
grue2_z = 0.05  # Grue Second Piece height
grue2_mass = 0.35  # Grue Second Piece mass
grue2_angle_value = [degrees_to_radian(150), degrees_to_radian(0)]
# Third Piece
grue3_x_value = [0.3, 0.4]  # Grue Third Piece length
grue3_z = 0.05  # Grue Third Piece height
grue3_mass = 0.6  # Grue Third Piece mass
grue3_angle_value = [degrees_to_radian(-60), degrees_to_radian(0)]

# -- Syringes --
#  todo: there are now negligees

# -- Windturbine --
windturbine_mass = 0.5  # Windturbine mass

# -- Counter-weight --
counterweight_x = 0.1
counterweight_position_x = - 0.1
counterweight_y = 0.1
counterweight_position_z = 0
counterweight_z = 0
counterweight_mass = 4

# -- Moving -- todo not now used
moving_max_x = 1
moving_max_z = 1

# -- Others --
I = 0.00001  # Inertia
g = 9.81  # [m / s²]
