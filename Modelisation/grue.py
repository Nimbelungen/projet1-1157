import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

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
# -- Mass of the system --
mass_sum = barge_mass + grue1_mass + grue2_mass + grue3_mass + grue4_mass + grapple_mass + windturbine_mass + \
           counterweight_mass

# -- Time --
step = 0.1  # dt [s]
end = 100.0  # [s]

# -- Numpy lists --
t = np.arange(0, end, step)  # [s] List of time
theta_rad = np.empty_like(t)  # [rad] List of Theta values
theta_deg = np.empty_like(t)  # [deg] List of theta value in degrees
omega_rad = np.empty_like(t)  # [rad/s] List of omega values
omega_deg = np.empty_like(t)  # [deg/s] List of omega values in degrees / seconds
E_k = np.empty_like(t)  # [J] List of kinetic energy values
E_g = np.empty_like(t)  # [J] List of gravitational energy

# -- Moving of the grue --
grue2_angle = np.empty_like(t)  # [rad] List of angles between grue piece 1 and grue piece 2
grue3_angle = np.empty_like(t)  # [rad] List of angles between grue piece 2 and grue piece 3
grue4_angle = np.empty_like(t)  # [rad] List of angles between grue piece 3 and grue piece 4
grue3_x = np.empty_like(t)  # [m] List of length of grue piece 3


def fill_array():
    """
    This function fills the lists of angles and lengths of the crane.
    :return: Nothing, it only modifies variables
    """
    grue2_angle[0] = grue2_angle_value[0]
    grue3_angle[0] = grue3_angle_value[0]
    grue4_angle[0] = grue4_angle_value[0]
    grue3_x[0] = grue3_x_value[0]
    step_grue2_angle = (grue2_angle_value[1] - grue2_angle_value[0]) / len(t)
    step_grue3_angle = (grue3_angle_value[1] - grue3_angle_value[0]) / len(t)
    step_grue4_angle = (grue4_angle_value[1] - grue4_angle_value[0]) / len(t)
    step_grue3_x = (grue3_x_value[1] - grue3_x_value[0]) / len(t)
    for i in range(len(t) - 1):
        grue2_angle[i + 1] = grue2_angle[i] + step_grue2_angle
        grue3_angle[i + 1] = grue3_angle[i] + step_grue3_angle
        grue4_angle[i + 1] = grue4_angle[i] + step_grue4_angle
        grue3_x[i + 1] = grue3_x[i] + step_grue3_x


# ---- Calculus Functions ---- Oder functions are in the 'formulas.py' file
def submersion_height():
    """
    Calculate the submerged height of the barge
    :return: If hc < barge height : the distance hc (for submerged height), where hc is the length follow the submerged
    z-axis of the barge. Otherwise, False
    """
    hc = mass_sum / (1000 * (barge_x * barge_y))
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
    return angle_max_x


def center_gravity(time):
    """
    This function calculates the coordinates of the center of gravity of the whole grue as a function of time.
    :type time: int
    :param time: Time. This is the index in the 'np' lists. These lists have been completed by the function fill_array()
    :return: A tuple that is the coordinate along the x and z axis of the center of gravity
    """
    # -- Barge --
    hc = submersion_height()
    hb = (barge_z / 2) - hc
    barge_cg = (0, hb)

    # -- Grue --
    # First Piece
    grue1_cg = (0, hb + (grue1_z / 2))

    # Second Piece
    grue2_cg_init = ((grue2_x / 2), (hb + (grue2_z / 2)))
    grue2_cg = rotate_coord(grue2_cg_init, grue2_angle[time])

    # Third Piece
    grue2_cg_init = (((math.cos(grue2_angle[time]) * grue2_x) + (grue3_x[time] / 2)),
                     (grue1_x + hb + (math.sin(grue2_angle[time]) * grue2_x) + (grue3_z / 2)))
    grue3_cg = rotate_coord(grue2_cg_init, grue3_angle[time])

    # Fourth Piece
    grue4_cg_init = (
        ((math.cos(grue2_angle[time]) * grue2_x) + (math.cos(grue3_angle[time]) * grue3_x[time]) + (grue4_x / 2)),
        (grue1_x + hb + (math.sin(grue2_angle[time]) * grue2_x) + (
                math.sin(grue3_angle[time]) * grue3_x[time]) + (grue4_z / 2)))
    grue4_cg = rotate_coord(grue4_cg_init, grue4_angle[time])

    # -- Syringes --
    # todo: there are now negligees

    # -- Windturbine -- = Coordonnées du bout de la partie 3 de la grue
    windturbine_cg_init = (
        ((math.cos(grue2_angle[time]) * grue2_x) + (math.cos(grue3_angle[time]) * grue3_x[time]) + grue4_x),
        (grue1_x + hb + (math.sin(grue2_angle[time]) * grue2_x) + (
                math.sin(grue3_angle[time]) * grue3_x[time]) + grue4_z))
    windturbine_cg = rotate_coord(windturbine_cg_init, grue4_angle[time])

    # -- Counterweight --
    counterweight_position_z = barge_z - hc
    counterweight_cg = (
        counterweight_position_x + (counterweight_x / 2), counterweight_position_z + (counterweight_y / 2))

    # Center of gravity
    cg = center_of_gravity_2d((barge_mass, barge_cg), (grue1_mass, grue1_cg), (grue2_mass, grue2_cg),
                              (grue3_mass, grue3_cg), (grue4_mass, grue4_cg), (windturbine_mass, windturbine_cg),
                              (counterweight_mass, counterweight_cg))
    return cg


def center_trust(angle):
    """
    Calculate the coordinate of the center of trust of the barge
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


def underwater_volume_mass(angle):  # todo: continue this function
    #  area of the trapeze times the width of the barge
    hc = submersion_height()
    parrallel_left = (hc - (math.tan(angle) * (barge_x / 2)))
    parrallel_right = (hc + (math.tan(angle) * (barge_x / 2)))
    area = ((parrallel_left * parrallel_right) * hc) / 2
    volume = area * barge_y

    return volume


def barge_inclination(time):
    """
    Binary search algorithm. It search the value of the angle of inclination of the barge.
    :type time: int
    :param time: Time. This is the index in the 'np' lists. These lists have been completed by the function fill_array()
    :return: The value of the angle in radians
    """
    global middle
    first = -math.pi / 2
    last = math.pi / 2  # todo: WARRING is the cause of some error
    find = False

    while first <= last and not find:
        middle = (first + last) / 2
        if abs(center_trust(middle)[0] - center_gravity(time)[0]) < 0.0000000000000001:
            # if center_trust(middle)[0] - center_gravity(time)[0] == 0:
            return middle
        else:
            if center_trust(middle)[0] < center_gravity(time)[0]:
                first = middle + 0.00000000000000001

            else:
                last = middle - 0.00000000000000001
    return middle


# ---- Simulation ----
def simulation():
    """
    This function completes the following lists (applying a binary search for angles)=
    - theta_rad and theta_deg \n
    - omega_rad and omega_deg
    - Energy : E_k, E_g and E_im
    :return: Nothing, it only modifies variables
    """
    # --- Theta lists ---
    # Rad
    for i in range(len(t)):
        theta_rad[i] = barge_inclination(i)
    # Deg
    for i_2 in range(len(t)):
        theta_deg[i_2] = rad_to_degrees(theta_rad[i_2])

    # --- Omega lists ---
    # Rad
    omega_rad[0] = None
    for j in range(len(t) - 1):
        omega_rad[j + 1] = (theta_rad[j + 1] - theta_rad[j]) / step
    # Deg
    for j_2 in range(len(t)):
        omega_deg[j_2] = rad_to_degrees(omega_rad[j_2])

    # --- Energy lists ---
    # Fill E_k List
    for k in range(len(t)):
        E_k[k] = I * ((omega_rad[k] ** 2) / 2)

    # Fill E_g List
    for m in range(len(t)):
        E_g[m] = mass_sum * g * center_gravity(m)[1]


def graph_angles():
    """
    This function creates the omega and theta graphs.
    :return: 2 graphs: the fist in radians and the second in degrees
    """
    # --- Max inclination value ---
    max_incl_rad = np.empty_like(t)
    max_incl_deg = np.empty_like(t)
    for i in range(len(t)):
        max_incl_rad[i] = maximum_inclination()
    for i_2 in range(len(t)):
        max_incl_deg[i_2] = rad_to_degrees(maximum_inclination())

    # --- Figure 1 : Radians ---
    plt.figure(1)
    plt.suptitle("Inclinaison et vitesse anglulaire [rad] et [rad/s]")
    plt.subplot(2, 1, 1)
    plt.plot(t, theta_rad, label="Thêta")
    plt.plot(t, max_incl_rad, label="Max angle", linestyle='dashed')
    plt.plot(t, -max_incl_rad, label="Max angle", linestyle='dashed')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(t, omega_rad, label="Omega")
    plt.legend()
    plt.show()

    # --- Figure 2 : Degrees ---
    plt.figure(2)
    plt.suptitle("Inclinaison et vitesse anglulaire [°] et [°/s]")
    plt.subplot(2, 1, 1)
    plt.plot(t, theta_deg, label="Thêta")
    plt.plot(t, max_incl_deg, label="Max angle", linestyle='dashed')
    plt.plot(t, -max_incl_deg, label="Max angle", linestyle='dashed')
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(t, omega_deg, label="Omega")
    plt.legend()
    plt.show()


def graph_energy():
    """
    This function creates the graph of energy
    :return: 2 graph, one with the 3 lines and another were the lines are separate
    """
    plt.figure(3)
    plt.suptitle("Energie [J]")
    plt.subplot(2, 1, 1)
    plt.plot(t, E_k, label="Énergie cinétique")
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(t, E_g, label="Énergie gravitationnel")
    plt.legend()
    plt.show()

    plt.figure(4)
    plt.suptitle("Energie [J]")
    plt.plot(t, E_k, label="Énergie cinétique")
    plt.plot(t, E_g, label="Énergie gravitationnel")
    plt.legend()
    plt.show()


# --- Lunch program ---
# -- Simulation and graphs
fill_array()
simulation()
graph_angles()
graph_energy()
# -- Print --
print("Simulation of the grue - group 11.57")
print(tabulate([["Information's about", "Radians", "Degrees"],
                ["Maximum Inclination", maximum_inclination(), rad_to_degrees(maximum_inclination())],
                ["Departure Inclination", barge_inclination(0), rad_to_degrees(barge_inclination(0))],
                ["Final Inclination", barge_inclination(-1), rad_to_degrees(barge_inclination(-1))]],
               headers="firstrow"))
print("Submersion Height = {}m".format(submersion_height()))
