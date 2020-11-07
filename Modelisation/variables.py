from formulas import *

# -- Barge --
barge_x = 0.5  # [m] Barge length
barge_y = 0.5  # [m] Barge width
barge_z = 0.3  # [m] Barge height
barge_mass = 0.3  # [kg] Barge mass

# -- Grue --

# - First Piece -
grue1_x = 0.15  # [m] Grue First Piece length
grue1_z = 0.05  # [m] Grue First Piece height
grue1_mass = 0.05  # [kg] Grue First Piece mass
# Grue First Piece angle is always 0

# - Second Piece -
grue2_x = 0.350  # [m] Grue Second Piece length
grue2_z = 0.05  # [m] Grue Second Piece height
grue2_mass = 0.35  # [kg] Grue Second Piece mass
grue2_angle_value = [degrees_to_radian(150), degrees_to_radian(0)]  # [deg] Angle of departure and arrival of this piece
# of grue. These angles are expressed as a function of the horizontal

# - Third Piece -
grue3_x_value = [0.3, 0.4]  # [m] Grue Third Piece start and finish length
grue3_z = 0.05  # [m] Grue Third Piece height
grue3_mass = 0.6  # [kg] Grue Third Piece mass
grue3_angle_value = [degrees_to_radian(-60), degrees_to_radian(0)]  # [deg] Angle of departure and arrival of this piece
# of grue. These angles are expressed as a function of the horizontal

# - Grapple -
grapple_mass = 0  # [kg] Mass of the grapple

# -- Syringes --
#  todo: there are now negligees

# -- Windturbine --
windturbine_mass = 0.5  # [kg] Windturbine mass
# The size of the turbine is not taken in consideration. It is considered to be in the grapple at all times.

# -- Counter-weight --
counterweight_x = 0.1  # [m] Counterweight length
counterweight_position_x = - 0.1  # [m] Position along the x-axis of the counterweight
counterweight_y = 0.1  # [m] Counterweight width
counterweight_z = 0  # [m] Counterweight height
counterweight_mass = 4  # [kg] Counterweight mass

# -- Moving -- todo not now used
moving_max_x = 0  # [m] Length of travel
moving_max_z = 0  # [m] Height of travel

# -- Others --
I = 0.00001  # [kg.m²]Inertia
g = 9.81  # [m / s²] Gravitational acceleration
