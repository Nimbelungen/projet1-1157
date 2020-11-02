import math
# You will find in this file the most important formulas that we will use to solve the different physics exercises.


# Center of gravity of a set of bodies
def center_of_gravity_3d(*args):
    """
    This function calculates the center of gravity of a set of n bodies.
    Each of these bodies has 2 coordinates (one x and one y). The system is thus in 2 dimensions.
    The formula that is use is : cg(x) = ( m1*d1(x) + m2*d2(x) + m3*d3(x) + ...) / (m1 + m2 + m3 + ...)
    (same form for the y axes).
    :type args: tuple
    :param args: Each argument is a tuple of values. They are of the form (mass, (x-coordinate, y-coordinate)).
    :return: A tuple witch the Coordinates in the form (x-coordinate, y-coordinate) of the center of gravity.

    NOT TESTED
    """
    mass_sum = 0
    for m in args:
        mass_sum += m[0]

    if mass_sum == 0:
        return False
    else:
        nominator_x = 0
        nominator_y = 0
        nominator_z = 0
        for x in args:
            mass_dist_x = x[0] * x[1][0]
            nominator_x += mass_dist_x
        for y in args:
            mass_dist_y = y[0] * y[1][1]
            nominator_y += mass_dist_y
        for z in args:
            mass_dist_z = z[0] * z[1][2]
            nominator_z += mass_dist_z

        cgx = nominator_x / mass_sum
        cgy = nominator_y / mass_sum
        cgz = nominator_z / mass_sum
        return tuple([cgx, cgy, cgz])


def rad_to_degrees(angle):
    """
    This function transforms a given angle in radians to angle in degrees.
    :type angle: float
    :param angle: The angle that we will 'change'
    :return: A float that is the angle in degrees
    """
    return (180 * angle) / math.pi


def degrees_to_radian(angle):
    """
    This function transforms a given angle in degrees to angle in radians.
    :type angle: float
    :param angle: The angle that we will 'change'
    :return: A float that is the angle in radian
    """
    return (angle * math.pi) / 180