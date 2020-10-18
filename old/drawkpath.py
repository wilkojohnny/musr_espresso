#!/usr/local/bin/python3.7

# drawkpath.py - draws a path between the points given. Meant to be for points in k-space but this is not a requirement.
# John Wilkinson 1/4/19

import numpy as np  # for numpy arrays
import argparse  # for parsing arguments from command line


def main():

    # get the out file location
    parser = argparse.ArgumentParser(description='Draw a path between a set of points')
    parser.add_argument("output_file_name", default="outfile.dat", help="name of the file to save the k-points to", nargs="?")
    arguments = parser.parse_args()

    # if there is an extension (just check for a dot!)
    outfile_location = arguments.output_file_name
    if "." not in outfile_location:
        outfile_location = outfile_location + ".dat"

    # get k points
    k_points = []
    while True:
        k, valid_k = get_k_vector()
        if valid_k:
            k_points.append(k)
        else:
            break

    # get how many points are between each pair of points
    n_points_line = 0
    while n_points_line == 0:
        print("How many points between each pair of points? > ", end="")
        try:
            n_points_line = abs(int(input()))
        except ValueError:
            n_points_line = 0

    k_path = []

    # open the output file
    outfile = open(outfile_location, 'w+')

    # for each consecutive pair in the array
    for i in range(0, len(k_points) - 1):
        # find the difference between this vector and the next one
        diff_k = k_points[i+1] - k_points[i]
        # divide this difference by the number of points per line
        dk = diff_k/n_points_line
        # and consecutively add this to k_points[i] and store in array and write to file
        k_j = k_points[i]
        for j in range(0, n_points_line):
            k_path.append(k_j)
            outfile.write(nparray_to_spacedstring(k_j) + " 0\n")  # the 0 "\n" is the weight for QE
            k_j = k_j + dk
    # don't forget the last point!
    k_path.append(k_points[-1])
    outfile.write(nparray_to_spacedstring(k_points[-1]) + " 0\n")

    # close the output file
    outfile.close()

    print(str(len(k_path)) + " points successfully written to " + outfile_location + ".")

    return 1


def get_k_vector() -> (np.ndarray, bool, ):
    print("k vector x,y,z > ", end="")
    try:
        # get the k vector input
        k_string = input()
        # try to parse the string - assume separated by commas
        k_str = k_string.split(",")
        k_x = float(k_str[0])
        k_y = float(k_str[1])
        k_z = float(k_str[2])
        k = np.array((k_x, k_y, k_z))
        return k, True
    except ValueError:
        return 0, False
    except IndexError:
        return 0, False


def nparray_to_spacedstring(inarray):
    # get a numpy array and convert it into a string of numbers with spaces in between
    outstring = ""
    for i in inarray:
        outstring = outstring + str(i) + " "
    # return the string (with a new line character for neatness)
    return outstring


if __name__ == '__main__':
    main()
