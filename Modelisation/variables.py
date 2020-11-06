# ---- Variables ----
# <name>_x = length [m]
# <name>_y = width  [m]
# <name>_z = height [m]
# <name>_mass = mass    [kg]

# -- Barge --
barge_x = 1  # Barge length
barge_y = 1  # Barge width
barge_z = 1  # Barge height
barge_mass = 1  # Barge mass

# -- Grue --
# First Piece
grue1_x = 1  # Grue First Piece length
grue1_z = 1  # Grue First Piece height
grue1_mass = 1  # Grue First Piece mass
# Grue First Piece angle is always 0
# Second Piece
grue2_x = 1  # Grue Second Piece length
grue2_z = 1  # Grue Second Piece height
grue2_mass = 1  # Grue Second Piece mass
grue2_angle_init = [0, 0]
# Third Piece
grue3_x = 1  # Grue Third Piece length
grue3_z = 1  # Grue Third Piece height
grue3_mass = 1  # Grue Third Piece mass

# -- Syringes --
syringes_mass = 1
syringes_d_grue1 = 1
syringes_d_grue21 = 1
syringes_d_grue22 = 1
syringes_d_grue3 = 1

# -- Windturbine --
windturbine_x = 1  # Windturbine length
windturbine_position_x = 1  # Windturbine distance center - axes system x => cg_x at time = 0
windturbine_z = 1  # Windturbine height
windturbine_mass = 1  # Windturbine mass

# -- Counter-weight --
counterweight_x = 0
counterweight_position_x = 0
counterweight_y = 0
counterweight_z = 0
counterweight_mass = 0

# -- Moving --
moving_max_x = 1
moving_max_z = 1

# -- Others --
I = 5  # Inertia
