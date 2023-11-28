def find_intersection(segment1, segment2):
    x1, y1 = segment1[0]
    x2, y2 = segment1[1]
    x3, y3 = segment2[0]
    x4, y4 = segment2[1]

    # Parametric equations for the lines
    denum = ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1))
    if denum != 0:
        ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denum
        ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denum

        # Check if the intersection point is within the parameter range of both segments
        if 0 < ua < 1 and 0 < ub < 1:
            # Calculate the intersection point
            intersection_x = x1 + ua * (x2 - x1)
            intersection_y = y1 + ua * (y2 - y1)
            return True, (intersection_x, intersection_y)
    
    # Segments do not intersect
    return False, None

def CollisionDetection(segments):
    """
    @ params:
    @ segments : a list of segments, each segment is defined as two points (its upper and lower bound)
    This function returns a (collision : boolean, data) : true if there is a collision between at least two of the rods, false otherwise
    """
    intersections = []

    for i in range(len(segments) - 1): # len(tige1)
        for j in range(i + 1, len(segments)): # len(tige2)
            collision, intersection_point = find_intersection(segments[i], segments[j])
            if collision:
                intersections.append(intersection_point)

    return (len(intersections)!=0), intersections
    """
    for seg1, seg2 in zip(tige1, tige2):
        collision, intersection_point = find_intersection(segments[i], segments[j])
            if collision:
                intersections.append(intersection_point)
    """




if __name__ == "__main__":
    # Example usage:
    segment1 = [(0, 0), (2, 2)]
    segment2 = [(1, 1), (3, 3)]
    segment3 = [(0, 2), (2, 0)]

    collision, intersection_point = CollisionDetection([segment1, segment3])

    if collision:
        print("Segments intersect at:", intersection_point)
    else:
        print("Segments do not intersect.")