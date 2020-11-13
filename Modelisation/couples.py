import matplotlib.pyplot as plt
import numpy as np

from Modelisation.variables import *
from formulas import *

# --- Simulation parameters
step = 0.001  # steps (dt) [s]
end = 20.0  # duration [s]
theta_0 = 0.0  # initial angle [rad]
omega_0 = 0.0  # initial angular velocity [rad/s]

t = np.arange(0, end, step)
theta = np.empty_like(t)
omega = np.empty_like(t)
a = np.empty_like(t)

# --- Moving of the grue ---
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


# --- Calculus Functions ---
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
                              (grue3_mass, grue3_cg), (grue4_mass, grue4_cg),
                              (windturbine_mass + grapple_mass, windturbine_cg),
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


def underwater_volume_mass(angle):
    """
    This function calculate the mass of the barge that is underwater
    :type angle: float
    :param angle: the angle of inclination of the barge
    :return: a float that is the mass of the underwater-barge. (a fraction of the total mass of the barge)
    """
    #  area of the trapeze times the width of the barge
    hc = submersion_height()
    parrallel_left = (hc - (math.tan(angle) * (barge_x / 2)))
    parrallel_right = (hc + (math.tan(angle) * (barge_x / 2)))
    area = ((parrallel_left + parrallel_right) * barge_z) / 2
    sub_volume = area * barge_y
    frac = sub_volume / (barge_x * barge_y * barge_z)
    mass_im = frac * barge_mass
    return mass_im


# --- Simulation and find angles ---
def simulation():
    # Initial conditions
    theta[0] = theta_0
    omega[0] = omega_0

    for i in range(len(t) - 1):
        dt = step

        # Calculation of torques
        C_g = -(force(mass_sum) * center_gravity(0)[0])  # We first consider the center of gravity in x as constant
        C_p = force(underwater_volume_mass(theta[i])) * center_trust(theta[i])[0]  # Depends on the angle
        C_d = -D * omega[i]
        C = C_g + C_p + C_d

        # Angle and angular velocity calculation
        a[i] = C / I
        omega[i + 1] = omega[i] + a[i] * dt  # todo: there is an error
        theta[i + 1] = theta[i] + omega[i + 1] * dt
        a[i + 1] = a[i]


def graphs():
    # Radians
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(t, theta)
    plt.plot([0, t[-1]], [maximum_inclination(), maximum_inclination()], '--r', label='Angle maximum')
    plt.plot([0, t[-1]], [-maximum_inclination(), -maximum_inclination()], '--r')
    plt.xlabel("t (s)")
    plt.ylabel("angle (rad)")
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(t, omega)
    plt.xlabel("t (s)")
    plt.ylabel("vistesse angulaire (rad/s)")
    plt.show()

    # Degrees
    plt.figure(2)
    plt.subplot(2, 1, 1)
    plt.plot(t, rad_to_degrees(theta))
    plt.plot([0, t[-1]], [rad_to_degrees(maximum_inclination()), rad_to_degrees(maximum_inclination())], '--r',
             label='Angle maximum')
    plt.plot([0, t[-1]], [-rad_to_degrees(maximum_inclination()), -rad_to_degrees(maximum_inclination())], '--r')
    plt.xlabel("t (s)")
    plt.ylabel("angle (°)")
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(t, rad_to_degrees(omega))
    plt.xlabel("t (s)")
    plt.ylabel("vistesse angulaire (°/s)")
    plt.show()

    # Phase diagram
    plt.figure(3)
    plt.title("Diagramme de phase")
    plt.plot(rad_to_degrees(omega), rad_to_degrees(theta))
    plt.xlabel("Theta (°)")
    plt.ylabel("Omega (°/s")
    plt.show()


# --- Lunch program ---
print(submersion_height())
fill_array()
simulation()
graphs()
