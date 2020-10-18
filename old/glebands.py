#!/usr/local/bin/python3.7

# glebands.py: Creates a GLE plot of the band structure which has been post-processed with bands.x
# John Wilkinson 29/3/19

import argparse  # for command line option parsing
import sys  # for yes/no option
import os  # for file size checking


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def ask_for_symmetry():
    # ask for the symmetry xvalues and their label: if no x input, then take as false
    sys.stdout.write("Symmetry x point > ")
    xposition_raw = input()

    # if the user didn't type in a number, end
    try:
        xpos = float(xposition_raw)
    except ValueError:
        return [False, False]

    # otherwise ask for a symmetry label
    sys.stdout.write("Symmetry label at x=" + xposition_raw + " > ")
    label = input()
    return [xpos, label]


def main():
    # get the .dat.gnu file from the input arguments (and while we're here, get all the others)
    parser = argparse.ArgumentParser(description='Plot a QE band structure with GLE')
    parser.add_argument("dat_gnu_file", help="location of the .dat.gnu file output by bands.x")
    parser.add_argument("--nogle", action="store_true", help="don't make the gle file")
    parser.add_argument("--noplot", action="store_true", help="don't plot the gle file")
    arguments = parser.parse_args()

    datgnufile_location = arguments.dat_gnu_file
    nogle = arguments.nogle
    noplot = arguments.noplot

    # do a quick check to see if the input file is sensible (i.e check it isn't empty, just has numbers, etc)
    # check the file exists, and that the extension is correct
    if datgnufile_location[-8:] != '.dat.gnu':
        print('File extension unexpected (wanted .dat.gnu)')
        if not query_yes_no('Proceed anyway?', "no"):
            print("Bye...")
            return
        else:
            datglefile_location = datgnufile_location + ".dat"
    else:
        datglefile_location = datgnufile_location[:-8] + "gle.dat"

    # try to open the file
    try:
        datgnufile = open(datgnufile_location, 'r')
    except IOError:
        print("File can't be opened. Bye...")
        return

    # now check the file has content
    if os.stat(datgnufile_location).st_size == 0:
        print('Input file has no content. Bye...')
        return

    # create a new file .dat, and put the data from the old file into it
    # create file
    try:
        datglefile = open(datglefile_location, 'w+')
    except IOError:
        print("Cannot create DAT file. Bye...")
        return

    # ask for the fermi level
    sys.stdout.write("Fermi/Highest occupied level (eV) > ")
    user_fermi = input()
    try:
        fermi = float(user_fermi)
    except ValueError:
        print("Fermi level taken to be 0 eV.")
        fermi = 0

    # make matrix for gle file like [[x1, band1, band2, band3],[x2, band1, band2, band3 ...
    datgle_matrix = []
    band_number = 1
    current_row = 0

    # for each line in the datgnufile
    try:
        datgnufile_lines = datgnufile.readlines()
        for datgnufile_line in datgnufile_lines:
            # get the values in the line
            values = list(map(float, datgnufile_line.split()))
            # check this line has something in it
            if len(values) > 0:
                xvalue = values[0]
                yvalue = values[1] - fermi
                # the bands.x code shows that the k values in the band will **always** be the same - so no need to check
                if band_number == 1:
                    # if this is the first band, do x
                    datgle_matrix.append([xvalue])
                # append on the band
                datgle_matrix[current_row].append(yvalue)

                # increment the row counter
                current_row += 1
            else:
                # if the line is blank, it means we're on a new band (if current row is >0)
                band_number += 1*(current_row > 0)
                current_row = 0
    except ValueError:
        # if we get here, there's a funny character somewhere in the file, so the format is probably wrong...
        print("File format not valid. Bye...")
        # close the opened files
        datgnufile.close()
        datglefile.close()
        return

    datgnufile.close()

    # save to a .dat file (datglefile) for input into gle
    datglefile.write("! k")
    for i in range(1, len(datgle_matrix[0])):
        datglefile.write(" Band" + str(i))

    for i in range(0, len(datgle_matrix)):
        datglefile.write("\n")
        for j in range(0, len(datgle_matrix[i])):
            datglefile.write(str(datgle_matrix[i][j]) + " ")

    datglefile.close()
    print("Written gle data to " + datglefile_location)

    if nogle:
        # if the user doesn't want the free GLE file, they can leave here
        return

    # for now - ask for symmetry points, and their labels (hopefully in the future use the bands.x symmetry)
    symmetry_points = []
    while True:
        currentsymmetrypoint = ask_for_symmetry()
        if currentsymmetrypoint[0] == currentsymmetrypoint[1] and currentsymmetrypoint[0] is False:
            break
        else:
            symmetry_points.append(currentsymmetrypoint)

    # make .gle file
    glefile_location = "bandplot.gle"
    # create file
    try:
        glefile = open(glefile_location, 'w+')
    except IOError:
        print("Cannot create GLE file. Exiting...")
        return

    # GLE preamble
    glefile.write("size 18 10\n\n")
    glefile.write("set font texcmr hei 0.3\n")
    glefile.write("begin graph\n")
    glefile.write("\tscale 0.8 0.75\n")
    glefile.write("\ttitle \"Band structure\"\n")
    glefile.write("\txaxis min 0 max " + str(datgle_matrix[-1][0]) + "\n")

    # do symmetry points in the x axis if any are defined
    if len(symmetry_points)>0:
        glefile.write("\txplaces ")
        for i in range(0, len(symmetry_points)):
            glefile.write(str(symmetry_points[i][0]) + " ")
        glefile.write("\n\txaxis nticks " + str(len(symmetry_points)))
        glefile.write("\n\txaxis grid")
        glefile.write("\n\txnames ")
        for i in range(0, len(symmetry_points)):
            glefile.write("\"" + str(symmetry_points[i][1]) + "\" ")
        glefile.write("\n")
        glefile.write("\txtitle \"k-symmetry point\"\n")
    else:
        glefile.write("\txtitle \"k-path position (arb. units)\"\n")

    glefile.write("\tytitle \"Energy (eV)\"\n")
    # tell GLE where the data to plot is
    glefile.write("\tdata \"" + datglefile_location + "\" ")
    for i in range(1, band_number):
        glefile.write("d" + str(i) + "=c1,c" + str(i+1) + " ")
    glefile.write("\n")
    for i in range(1, band_number):
        glefile.write("\t d" + str(i) + " line \n")
    glefile.write("end graph")

    glefile.close()
    print("Written GLE file to " + glefile_location)

    # unless asked not to with --noplot, plot the .gle file with gle
    if not noplot:
        os.system("gle " + glefile_location)


if __name__ == "__main__":
    main()
