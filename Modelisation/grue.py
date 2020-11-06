import matplotlib.pyplot as plt
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
mass_sum = windturbine_mass + barge_mass + grue1_mass + grue2_mass + grue3_mass + counterweight_mass

# Time
step = 1  # dt [s]
end = 60.0  # [s]

# -- Numpy lists --
t = np.arange(0, end, step)  # List of time
theta = np.empty_like(t)
omega = np.empty_like(t)
E_k = np.empty_like(t)

# WARNING
grue2_angle = np.empty_like(t)
grue3_angle = np.empty_like(t)
grue3_x = np.empty_like(t)  # todo: attention nom de la variable


def fill_array():
    grue2_angle[0] = grue2_angle_value[0]
    grue3_angle[0] = grue3_angle_value[0]
    grue3_x[0] = grue3_x_value[0]
    step_grue2_angle = (grue2_angle_value[1] - grue2_angle_value[0]) / len(t)
    step_grue3_angle = (grue3_angle_value[1] - grue3_angle_value[0]) / len(t)
    step_grue3_x = (grue3_x_value[1] - grue3_x_value[0]) / len(t)
    for i in range(1, len(t)):
        grue2_angle[i] = grue2_angle[i - 1] + step_grue2_angle
        grue3_angle[i] = grue3_angle[i - 1] + step_grue3_angle
        grue3_x[i] = grue3_x[i - 1] + step_grue3_x


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


def center_gravity(time):
    # -- Barge --
    hc = submersion_height()
    hb = (barge_z / 2) - hc
    barge_cg = (0, hb)

    # -- Grue --
    # First Piece OK
    grue1_cg = (0, hb + (grue1_z / 2))

    # Second Piece OK
    grue2_cg_init = ((grue2_x / 2), (hb + (grue2_z / 2)))
    grue2_cg = rotate_coord(grue2_cg_init, grue2_angle[time])

    # Third Piece OK
    grue2_cg_init = (((math.cos(grue2_angle[time]) * grue2_x) + (grue3_x[time] / 2)),
                     ((math.sin(grue2_angle[time]) * grue2_x) + (grue3_z[time] / 2)))
    grue3_cg = rotate_coord(grue2_cg_init, grue3_angle[time])

    # -- Syringes --
    # todo: there are now negligees

    # -- Windturbine -- = Coordonnées du bout de la partie 3 de la grue
    windturbine_cg_init = (((math.cos(grue2_angle[time]) * grue2_x) + (grue3_x[time])),
                           ((math.sin(grue2_angle[time]) * grue2_x) + (grue3_z[time])))
    windturbine_cg = rotate_coord(windturbine_cg_init, grue3_angle[time])

    # -- Counterweight --
    counterweight_cg = (counterweight_position_x, (counterweight_y / 2))

    # Center of gravity
    cg = center_of_gravity_2d((barge_mass, barge_cg), (grue1_mass, grue1_cg), (grue2_mass, grue2_cg),
                              (grue3_mass, grue3_cg), (windturbine_mass, windturbine_cg),
                              (counterweight_mass, counterweight_cg))
    return cg


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


def barge_inclination(time):
    """
    Binary search algorithm. It search the value of the angle of inclination of the barge.
    :return: The value of the angle in radians
    """
    first = 0.0
    last = maximum_inclination()
    find = False

    while first <= last and not find:
        middle = (first + last) / 2
        if center_trust(middle)[0] == center_gravity(time)[0]:
            return middle
        else:
            if center_trust(middle)[0] < center_gravity(time)[0]:
                first = middle + 0.00000000000000001

            else:
                last = middle - 0.00000000000000001


# ---- Simulation ----
def simulation():
    # Fill theta list
    for i in range(len(t)):
        theta[i] = barge_inclination(t)
    omega[0] = 0
    # Fill omega list
    for j in range(len(t) - 1):
        omega[j + 1] = (theta[j + 1] - theta[j]) / step
    # Fill E_k List
    for k in range(len(t)):
        E_k[k] = I * ((omega[k] ** 2) / 2)


def graphique_angles():
    plt.figure(1)

    plt.subplot(2, 1, 1)
    plt.plot(t, theta, label="Thêta")
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(t, omega, label="Omega")
    plt.legend()

    plt.show()
