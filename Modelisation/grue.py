import math

from Modelisation.variables import *

""" COORDINATE SYSTEM
The following code takes place in a 3-dimensional coordinate system. 
The X axis is horizontal (length)
The Y-axis is horizontal (width)
The Z axis is vertical (height)
The origin is positioned in the middle of the barge along the X and Y axis and at water level along the Z axis.
"""


# ---- Calculus Functions ----
def submersion_height():
    """
    Calculate the submerged height of the barge
    :return: If hc < barge height : the distance hc (for submerged
    height), where hc is the length follow the submerged z-axis of the barge. Otherwise, False
    """
    sub_volume = (
                             mass_mass + barge_mass + grue1_mass + grue2_mass + grue3_mass + syringe_12_mass + syringe_23_mass) / 1000
    hc = sub_volume / (barge_x * barge_y)
    if hc < barge_z:
        return hc
    else:
        return False


def maximum_inclination():
    """
    This function calculates the maximum tilt angles before the barge sinks along the X- and Y-axis.
    :return: A tuple. It contains the value in radians of the angle first along the X-axis and then along the Y-axis.
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

    # --- Along th Y axis ---
    # Fist method
    tan_theta1_y = (barge_z - submersion_height()) / (barge_y / 2)
    angle_max1_y = math.atan(tan_theta1_y)
    # Second method
    tan_theta2_y = submersion_height() / (barge_y / 2)
    angle_max2_y = math.atan(tan_theta2_y)
    # Angle
    if angle_max1_y <= angle_max2_y:
        angle_max_y = angle_max1_y
    else:
        angle_max_y = angle_max2_y

    # --- Return ---
    return tuple([angle_max_x, angle_max_y])


def center_gravity(init):
    hc = submersion_height()
    dist_axe_x_center_gravity_barge = (barge_y / 2) - hc
    hd = barge_y - hc

    # Center gravity Barge
    barge_cg = (0, 0, dist_axe_x_center_gravity_barge)

    # Center gravity Grue
    grue1_cg = (0, 0, (grue1_y / 2) + hd)
    if init:
        grue2_cg = (0, 0, "z")
        grue3_cg = (0, 0, "z")
    else:
        grue2_cg = (0, 0, "z")
        grue3_cg = (0, 0, "z")

    # Center gravity Mass
    if init:
        "position de la mass à prendre"
    else:
        "position de la masse à déposer"

    # Center gravity Syringes todo: in function of the grue ?
    "position de la première seringue ne change pas"
    if init:
        "position de la deuxième seringe"
    else:
        "position de la deuxième seringue après déplacement"
