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


# ---- Calculus Functions ----
# Oder functions are in the 'formulas.py' file
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


def center_gravity(init):
    """
    Calculate the center of gravity of a set of body's
    :type init: bool
    :param init: True if it's the initial situation. Otherwise, False
    :return: A tuple with the coordinate along X- and Z-axis of the center of gravity of all the variable <name>_mass
    """
    hc = submersion_height()
    dist_axe_x_center_gravity_barge = (barge_y / 2) - hc
    hd = barge_y - hc

    # Center gravity Barge
    barge_cg = (0, dist_axe_x_center_gravity_barge)

    # Center gravity Grue
    grue1_cg = (0, (grue1_z / 2) + hd)
    if init:
        grue2_cg = (0, "z")
        grue3_cg = (0, "z")
    else:
        grue2_cg = (0, "z")
        grue3_cg = (0, "z")

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

    return tuple(cgx, cgz)


def center_trust(init, angle):
    """
    Calculate the coordinate of the center of trust of a trapezium if init == False
    :type init: bool
    :type angle: float
    :param init: True if it's the initial situation. Otherwise, False
    :param angle: The angle of inclination that the barge undergoes, changing the coordinate system and causing the
    submerged shape change.
    :return: A tuple with the coordinate along X- and Z-axis of the center of trust
    """
    hc = submersion_height()

    if init:
        ctx = barge_x / 2
        cty = barge_y / 2
        ctz = hc / 2
        return tuple([ctx, cty, ctz])
    else:
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


# --- Find the inclination of the barge ---
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
