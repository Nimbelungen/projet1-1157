import numpy as np

from Modelisation.variables import *
from formulas import *

"""COORDINATE SYSTEM 
The following code takes place in a 3-dimensional coordinate system. However, some dimensions will be regularly ignored
(especially the Y component). A tuple with 2 coordinates is thus composed of the x-y coordinates. 
The X axis is horizontal (length) 
The Y-axis is horizontal (width) T
he Z axis is vertical (height) 
The origin is positioned in the middle of the barge along the X and Y axis and at water level along the Z axis. 
"""

# ---- Calculus ----
# Mass of the system
mass_sum = windturbine_mass + barge_mass + grue1_mass + grue2_mass + grue3_mass + (2 * syringes_mass) + \
           counterweight_mass

# Time
step = 1  # dt [s]
end = 50.0  # [s]

# Numpy list
t = np.arange(0, end, step)  # List of time
d_x = np.arange(0, moving_max_x, (moving_max_x / len(t)))  # List of x displacement
d_y = np.arange(0, moving_max_z, (moving_max_z / len(t)))  # List of x displacement

# Syringes
alpha = np.empty_like(t)  # Angle between pieces 1 and 2 of the grue
beta = np.empty_like(t)  # Angle between pieces 2 and 3 of te barge
syringes_length = np.arange(90, 155, (65 / len(t)))

# Grue 2 length
grue2_length = np.empty_like(t)


# ---- Simulation Functions ----
def angle_syringes_grue():
    # todo: link the time and the length of the syringes
    # todo: link the angles of the grue parts and the length of the syringes
    # Necessary condition : moving_max_x <= grue2_length[t max] + grue3_x
    if moving_max_x <= grue2_length[-1] + grue3_x:
        pass
    else:
        print("REALITY ERROR: The Grue is too short")
        raise


# ---- Calculus Functions ---- Oder functions are in the 'formulas.py' file
def submersion_height():
    """
    Calculate the submerged height of the barge
    :return: If hc < barge height : the distance hc (for submerged
    height), where hc is the length follow the submerged z-axis of the barge. Otherwise, False
    """
    sub_volume = mass_sum / 1000
    hc = sub_volume / (barge_x * barge_y)
    if hc < barge_z:
        return hc
    else:
        return False


def maximum_inclination():
    """
    This function calculates the maximum tilt angles before the barge sinks along the X-axis.
    :return: The value in radians of the angle along the X-axis.
    """
    # --- Along the X-axis ---
    # Fist method
    tan_theta1_x = (barge_z - submersion_height()) / (barge_x / 2)
    angle_max1_x = math.atan(tan_theta1_x)
    # Second method
    tan_theta2_x = submersion_height() / (barge_x / 2)
    angle_max2_x = math.atan(tan_theta2_x)
    # Angle
    if angle_max1_x <= angle_max2_x:
        angle_max_x = angle_max1_x
    else:
        angle_max_x = angle_max2_x

    # --- Return ---
    return tuple([angle_max_x])


def center_gravity(t):
    # -- Barge --
    hc = submersion_height()
    hb = (barge_z / 2) - hc
    barge_cg = (0, hb)

    # -- Grue --
    # First Piece
    grue1_cg = (0, hb + (grue1_z / 2))
    # Second Piece
    grue2_cg_init = ((grue2_x / 2), (grue2_z / 2))
    # todo: with the angle that produce the syringe
    # Third Piece

    # -- Syringes --

    # -- Windturbine --
    windturbine_cg = (windturbine_position_x, (windturbine_z / 2))  # todo add the x and y displacement

    # -- Counterweight --
    counterweight_cg = (counterweight_position_x, (counterweight_y / 2))


def center_trust(angle):
    """
    Calculate the coordinate of the center of trust of a trapezium if init == False
    :type angle: float
    :param angle: The angle of inclination that the barge undergoes, changing the coordinate system and causing the
    submerged shape change.
    :return: A tuple with the coordinate along X- and Z-axis of the center of trust
    """
    hc = submersion_height()

    # --- The coordinate system also rotates ---
    """ Formulas
    lxc = lctx = ( l * (h1 + 2 * h2) ) / (3 * (h1 + h2) )
    hc = lctz = (h1 ** 2 + h1 * h2 + h2 ** 2) / (3 * (h1 + h2) )
    Where :
     -  l = barge_x
     - h1 = parrallel_right 
     - h2 = parrallel_left 
       => the trapeze is in the wrong direction in the slides
    """
    # Length of the 2 parallel sides of the trapeze
    parrallel_left = (hc - (math.tan(angle) * (barge_x / 2)))
    parrallel_right = (hc + (math.tan(angle) * (barge_x / 2)))

    # Application of formulas
    lctx = (barge_x / 3) * ((parrallel_right + (2 * parrallel_left)) / (parrallel_right + parrallel_left))
    lctz = ((parrallel_right ** 2) + (parrallel_right * parrallel_left) + (parrallel_left ** 2)) / (
            3 * (parrallel_right + parrallel_left))

    # Coordinates in the inclined system
    ctx_rotate = (barge_x / 2) - lctx
    ctz_rotate = hc - lctz

    # --- Rotation of the system ---
    ctx = (ctx_rotate * math.cos(-angle)) - (ctz_rotate * math.sin(-angle))
    ctz = (ctx_rotate * math.sin(-angle)) + (ctz_rotate * math.cos(-angle))

    return tuple([ctx, ctz])


def barge_inclination():
    """
    Binary search algorithm. It search the value of the angle of inclination of the barge.
    :return: The value of the angle in radians
    """
    first = 0.0
    last = maximum_inclination()
    find = False

    while first <= last and not find:
        middle = (first + last) / 2
        if center_trust(False, middle)[0] == center_gravity(False)[0]:
            return middle
        else:
            if center_trust(False, middle)[0] < center_gravity(False)[0]:
                first = middle + 0.00000000000000001

            else:
                last = middle - 0.00000000000000001

