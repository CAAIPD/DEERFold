import numpy as np

def direction_vector(P1, P2):
    """Calculate the direction vector between two points"""
    return P2 - P1

def angle_between_vectors(v1, v2):
    """Calculate the angle between two vectors (in radians)"""
    cos_theta = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)
    return np.arccos(cos_theta)

def find_min_angle_pair(A, B, C, D):
    """
    Calculate the angles between two pairs of line segments (AD & BC, AC & BD)
    and return the pair with the smallest angle and the angle value.

    Parameters:
        A, B, C, D: 3D points (numpy arrays)

    Returns:
        min_pair: The pair of segments with the smallest angle
        min_angle_rad: The smallest angle (in radians)
        min_angle_deg: The smallest angle (in degrees)
    """
    # Calculate direction vectors for each segment
    AD = direction_vector(A, D)
    BC = direction_vector(B, C)
    AC = direction_vector(A, C)
    BD = direction_vector(B, D)
    
    # Calculate angles between segment pairs
    angles = {
        '1,4,2,3': angle_between_vectors(AD, BC),
        '1,3,2,4': angle_between_vectors(AC, BD)
    }
    
    # Find the pair with the smallest angle
    min_pair = min(angles, key=angles.get)
    min_angle_rad = angles[min_pair]
    min_angle_deg = np.degrees(min_angle_rad)
    
    return min_pair

# A = np.array([5.888,   75.570,  46.012])
# B = np.array([-4.476,  56.428,  35.286])
# C = np.array([-11.041, 58.550,  27.901])
# D = np.array([-7.576,  56.425,  24.457])
# min_pair = find_min_angle_pair(A, B, C, D)
# print(min_pair)