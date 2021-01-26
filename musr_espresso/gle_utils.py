# gle_utils.py - utilities for plotting with GLE
import subprocess  # allows us to run shell commands

color_list = ["black", "red", "green", "blue", "teal", "gold", "magenta", "mediumpurple"]
marker = ["wcircle", "wtriangle", "wsquare", "wdiamond", "odot", "cross", "flower", "club", "heart", "spade", "phone"]


# plot simple x y data
def plot_xy(xdata, ydata, title='xy data', xtitle="xtitle", ytitle="ytitle", xrange=None, yrange=None, size=None,
            postamble="", legend=None, colors=None):
    # fix mutable and null arguments
    if size is None:
        size = [10, 10]
    if legend is None:
        legend = "data"
    if colors is None:
        colors = color_list
    # make the inputs lists of lists if not already
    if not isinstance(xdata[0], list):
        xdata = [xdata]
    if not isinstance(ydata[0], list):
        ydata = [ydata]
    # save the x y data into a text file
    data_file = open("plot.dat", "w")

    # get size of legend
    if isinstance(legend, list):
        legend_length = len(legend)
    else:
        legend_length = 1

    # write the title (based on the legend)
    data_file.write("! ")
    if legend is not None and isinstance(legend, list):
        for legend_entry in legend:
            data_file.write(legend_entry + "_x " + legend_entry + "_y ")
    else:
        data_file.write("x y")
    data_file.write("\n")

    for ix in range(0, len(xdata[0])):
        data_file_line = ""
        for i_series in range(0, len(xdata)):
            try:
                data_file_line = data_file_line + str(xdata[i_series][ix]) + " " + str(ydata[i_series][ix]) + " "
            except IndexError:
                # if we are out of range of an array (i.e due to different sized arrays), just write a *
                data_file_line = data_file_line + '* * '
        data_file.write(data_file_line + "\n")

    data_file.close()

    # write the gle file
    gle_file = open("plot.gle", "w")

    gle_file.write("size " + str(size[0]) + " " + str(size[1]) + "\n")
    gle_file.write("set texlabels 1\n")
    gle_file.write("begin graph \n \t size " + str(size[0]) + " " + str(size[1]) + "\n")

    data_line = "\t "
    linestyle_lines = []
    # for each series on the graph
    for i_series in range(1, legend_length+1):
        if legend_length == 1:
            legend_entry = legend
        else:
            legend_entry = legend[i_series-1]
        # tell gle which columns this is
        data_line = data_line + "d" + str(i_series) + "=c" + str(2*i_series-1) + ",c" + str(2*i_series) + " "
        # style the series
        linestyle_lines.append("\t d" + str(i_series) + " line marker fcircle msize 0.11 color "
                               + color_list[i_series-1] + " key \"" + legend_entry + "\"\n")

    gle_file.write("\t data \"plot.dat\" " + data_line + " \n")
    gle_file.writelines(linestyle_lines)
    gle_file.write("\t title \"" + title + "\" \n\t xtitle \"" + xtitle + "\" \n\t ytitle \"" + ytitle + "\" \n")

    # axis limits
    if xrange is not None:
        gle_file.write("\t xaxis ")
        if xrange[0] is not None:
            gle_file.write("min " + str(xrange[0]) + " ")
        if xrange[1] is not None:
            gle_file.write("max " + str(xrange[1]) + " ")
        gle_file.write("\n")

    if yrange is not None:
        gle_file.write("\t yaxis ")
        if yrange[0] is not None:
            gle_file.write("min " + str(yrange[0]) + " ")
        if xrange[1] is not None:
            gle_file.write("max " + str(yrange[1]) + " ")
        gle_file.write("\n")

    gle_file.write("end graph \n")

    gle_file.write(postamble)

    gle_file.close()

    # now run the gle command
    gle_output = subprocess.run("gle plot.gle", shell=True, check=True)
    if gle_output.returncode == 0:
        # if successful finish, show the plot
        subprocess.Popen(["gv", "plot.eps"], stdout=subprocess.PIPE)
    else:
        print("error running gle.")


def plot_scatter(data_files, title='xydata', legend=None, file_location='plot.gle'):
    # find out how many data_files there are
    if isinstance(data_files, str):
        # only one data file
        nData = 1
    else:
        # multiple data files
        nData = len(data_files)

    gle_file = open(file_location, "w")

    gle_file.write("size 10 10\n\n")
    gle_file.write("amove 1.4 1\n")
    gle_file.write("begin graph\n")
    gle_file.write("\tsize 8 8\n")
    gle_file.write("\tfullsize\n\n")

    if nData == 1:
        gle_file.write("\tdata " + data_files + "\n")
        gle_file.write("\td1 marker " + marker[0] + " msize 0.2 color blue\n")
    else:
        for i in range(0, nData):
            gle_file.write("\tdata " + data_files[i] + "\n")
            gle_file.write("\td" + str(i+1) + " marker " + marker[i] + " msize 0.2 color " + color_list[i]
                           + " key " + legend[i] + "\n")

    # draw dotted 0 line
    gle_file.write("\tlet d" + str(nData+1) + " = 0 from 0 to 10\n")
    gle_file.write("\td" + str(nData+1) + " line lstyle 2\n")

    gle_file.write("\n")
    gle_file.write("\ttitle \"" + title + "\"\n")

    gle_file.write("end graph\n")

    gle_file.close()


def plot_bands(dat_gnu_file_location: str, symm_x_labels: list = None, fermi_energy: float = 0):
    """
    plot a band structure generated with bands.x with GLE
    :param dat_gnu_file_location: location of the gnuplot data file to convert into gle-ish
    :param symm_x_labels: array of [[symm_x position, symm_label]] to plot
    :param fermi_energy: energy of the fermi / higest occupied level
    :return:
    """

    # make matrix for gle file like [[x1, band1, band2, band3],[x2, band1, band2, band3 ...
    datgle_matrix = []
    band_number = 1
    current_row = 0

    with open(dat_gnu_file_location, 'r') as datgnufile:
        # for each line in the datgnufile
        try:
            datgnufile_lines = datgnufile.readlines()
            for datgnufile_line in datgnufile_lines:
                # get the values in the line
                values = list(map(float, datgnufile_line.split()))
                # check this line has something in it
                if len(values) > 0:
                    xvalue = values[0]
                    yvalue = values[1] - fermi_energy
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
            return

    # sort out the file names
    datglefile_location = dat_gnu_file_location[0:-4]
    glefile_location = ''
    if not datglefile_location.endswith('.dat'):
        glefile_location += '.gle'
        datglefile_location += '.dat'
    else:
        glefile_location = datglefile_location[0:-4] + '.gle'

    with open(datglefile_location, 'w') as datglefile:
        # save to a .dat file (datglefile) for input into gle
        datglefile.write("! k")
        for i in range(1, len(datgle_matrix[0])):
            datglefile.write(" Band" + str(i))

        for i in range(0, len(datgle_matrix)):
            datglefile.write("\n")
            for j in range(0, len(datgle_matrix[i])):
                datglefile.write(str(datgle_matrix[i][j]) + " ")

    print("Written gle data to " + datglefile_location)

    # make .gle file
    with open(glefile_location, 'w+') as glefile:
        # GLE preamble
        glefile.write("size 18 10\n\n")
        glefile.write("set font texcmr hei 0.3\n")
        glefile.write("begin graph\n")
        glefile.write("\tscale 0.8 0.75\n")
        glefile.write("\ttitle \"Band structure\"\n")
        glefile.write("\txaxis min 0 max " + str(datgle_matrix[-1][0]) + "\n")

        # do symmetry points in the x axis if any are defined
        if len(symm_x_labels)>0:
            glefile.write("\txplaces ")
            for i in range(0, len(symm_x_labels)):
                glefile.write(str(symm_x_labels[i][0]) + " ")
            glefile.write("\n\txaxis nticks " + str(len(symm_x_labels)))
            glefile.write("\n\txaxis grid")
            glefile.write("\n\txnames ")
            for i in range(0, len(symm_x_labels)):
                glefile.write("\"" + str(symm_x_labels[i][1]) + "\" ")
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

    print("Written GLE file to " + glefile_location)

    # unless asked not to with --noplot, plot the .gle file with gle
    subprocess.run("gle " + glefile_location, shell=True)

    return True
