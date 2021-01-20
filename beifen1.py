'''
[Intro]
Plot
# (1) the overlap of each xi window,
# (2) the variance in each trajectory,
# (3) potential of mean force,
# (4) recrossing factor,
for a single task.

[Usage]
run `python <this file name>` then you will be ask to input a path containing a result for a single task.
Then you'll get all four figures above if your task ended normally.

Attention, please! The former figures will be deleted when the program started running.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as itp # Interpolate
import os

color = ['#00447c', '#ae0d16', '#47872c', '#800964']
# SHU Blue, Weichang Red, Willow Green, SHU Purple
# This color scheme can be easily obtained on the official website `vi.shu.edu.cn`.

def input_path():
    path = input('Please input the result folder: ')
    return path

def plot_parameters(title):
    print('Plotting '+title)
    plt.figure(figsize=(4, 3))
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix" # The font closet to Times New Roman
    # plt.rcParams["mathtext.fontset"] = "stix" # The font closet to Times New Roman
    plt.rcParams['xtick.direction'] = 'in'

    file = title+'.png'
    if os.path.exists(file):
        os.remove(file)
def plot_save(name):
    plt.tight_layout()
    plt.savefig(name + '.png', dpi=600) # .svg is recommended!
    plt.clf()
    plt.close()

def get_tail(fname):
    fp = open(fname, 'r')
    last_line = fp.readlines()[-1]
    fp.close()
    return last_line
def get_xav(fname):
    line = get_tail(fname)
    line_element = line.split()
    Nl = len(line_element)
    if Nl < 3:
        raise ValueError("The last line doesn't contain 'mean' and 'xav2'. ")
    elif line_element[3] == "=" or line_element[3] == "*":
        print("error\n",fname)
    else:
        mean = np.float(line_element[3])
        xav2 = np.float(line_element[4])
    return mean, xav2
def my_gaussian(x, mean, xav2):
    y=(1.0/np.sqrt(2.0*np.pi*xav2))*np.exp(-(x-mean)**2/(2.0*xav2))
    return y

def get_xilist(path):
    fileList = os.listdir(path)
    xiList = []
    for file in fileList:
        if file[0:17] == 'umbrella_sampling':
            xiList.append(np.float(file[18:25].rstrip('.')))
    xiList.sort()
    return xiList

def plot_overlap(path, xiList):
    a=0
    title = 'Overlap'
    plot_parameters(title)

    # plt.yticks([])# No ticks and labels in y axis

    resolution = 2000
    extend = 0 # 3E-2

    xiMin = np.min(xiList)
    xiMax = np.max(xiList)
    length = len(xiList)


    x_new = np.linspace(xiMin - extend, xiMax + extend, resolution)
    y_sum = np.zeros(resolution) # Total density line
    A=[]
    for i in range(length):
        if a>=4:
            a=0
        fname = path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[i])
        if os.path.isfile(fname):
            # Gaussian smearing
            mean, xav2 = get_xav(fname)
            y_new = my_gaussian(x_new, mean, xav2)

            # Find biggest population
            maxPop = max(y_new)
            A.append(maxPop)
            maxA=max(A)

            y_sum += y_new
            # sum all population
            if xav2 > 5.0E-5:
                print("[WARNING] May be too various in xi = {0:.4f}! ".format(xiList[i]))
                plt.plot(x_new, y_new, lw=2, c=color[1])  #alpha=0.8)   #red
            else:
                plt.plot(x_new, y_new, lw=0.5, c=color[a],  alpha=.6)
                a+=1

    # Plot summation and difference # blue
    # plt.plot(x_new, y_sum, lw=1, c=color[0], label='Summation of all populations') # SHU Blue

    # plt.xlabel('Reaction Coordinate / Å')
    plt.xlabel('Reaction Coordinate')
    plt.ylabel('Population')
    # plt.legend(loc="best") # No legend

    plt.xlim(xiMin - extend, xiMax + extend)
    plt.ylim(0, maxA*1.1)
    # plt.show()

    plot_save(title)
def plot_variance(path, xiList):
    title = 'Variance'
    plot_parameters(title)

    xiMin = np.min(xiList)
    xiMax = np.max(xiList)
    length = len(xiList)

    for i in range(length):
        fname = path + "/umbrella_sampling_{0:.4f}.dat".format(xiList[i])
        f = open(fname, 'r')
        fLines = f.readlines()
        #print ("a00",fLines)
        f.close()

        timeEvolution = []
        xivar = []
        for line in fLines[15:]:
            lineSplit = line.split()
            timeEvolution.append(np.float(lineSplit[2]))
            xivar.append(np.float(lineSplit[-1]))

        # timeStep = np.float(fLines[6].split()[3])
        timeEvolution = [x * 1E-6 for x in timeEvolution]
        # xivarDelta = []
        # for i in range(len(xivar)-1):
        #     xivarDelta.append(np.abs(xivar[i+1] - xivar[i]))
        # x = range(len(xivar))
        # plt.yscale('log')

        # # Shifted
        # for i in range(len(xivar)):
        #     xivar[i] = xivar[i] - xivar[0]

        plt.xlabel('$t$ / ns')
        plt.ylabel('Variance')

        # Scientific notation for y axis
        # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        from matplotlib import ticker
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        plt.gca().yaxis.set_major_formatter(formatter)

        # color = (np.random.rand(), np.random.rand(), np.random.rand())
        xivar1=np.array(xivar)
        A=np.std(xivar1)
        if A>=1e-5:
            print (fname)
        plt.plot(timeEvolution, xivar, lw=0.2, c=color[0], alpha=0.5)

    plot_save(title)
def plot_pmf(path):
    title = 'PMF'
    plot_parameters(title)

    try:
        f = open(path + '/potential_of_mean_force.dat', 'r')
    except FileNotFoundError:
        print('[ERROR] {} file not found! '.format(title))
    else:
        fLines = f.readlines()
        f.close()

        xi = []
        pmf = []
        for i in fLines[12:-1]:
            xi.append(np.float(i.split()[0]))
            pmf.append(np.float(i.split()[1]))
        N=len(pmf)
        pmf=pmf-np.ones(N)*pmf[0]
        print (max(pmf))

        # Let W(xi=0) = 0!
        # xiAbs = np.abs(xi)
        # xiZeroIndex = list(xiAbs).index(min(xiAbs))
        # pmf = [x - pmf[xiZeroIndex] for x in pmf]

        plt.xlabel(r'$\xi$')
        plt.ylabel(r'$W(\xi)$ / eV')

        # Choose the fit range of this plot
        # plt.xlim(xi[0], xi[-1])
        # plt.ylim(pmf[0]-max(pmf)*0.1, max(pmf)*1.1) #The upper adds 0.1*max(pmf), the lower region should be added the same.

        plt.plot(xi, pmf, c=color[0])
        plt.tick_params(direction='in')
        # # plot a zoomed subfigure
        # xiMaxIndex = pmf.index(max(pmf)) # the position of maximum
        # extend = 80 # Extra points to plot
        # ximax = xi[xiMaxIndex - extend : xiMaxIndex + extend]
        # pmfmax = pmf[xiMaxIndex - extend : xiMaxIndex + extend]
        #
        # # # Find maximum of xi
        # f = itp(ximax, pmfmax, k=4)
        # cr_pts = f.derivative().roots()
        # cr_pts = np.append(cr_pts, (ximax[0], ximax[-1]))  # also check the endpoints of the interval
        # cr_vals = f(cr_pts)
        # #! min_index = np.argmin(cr_vals)
        # max_index = np.argmax(cr_vals)
        # pmfMax, xiMax = cr_vals[max_index], cr_pts[max_index]
        #
        # subfig = plt.axes([.3, .5, .5, .4])
        # subfig.plot(ximax, pmfmax, c=color[0])
        #
        # subfig.axvline(x=xiMax, c=color[0], lw=0.5, linestyle='--')
        # subfig.axhline(y=pmfMax, c=color[0], lw=0.5, linestyle='--')
        #
        # plt.setp(subfig, xlim=[min(ximax), max(ximax)])
        plt.show()
        plot_save(title)
def plot_rexFactor(path):
    title = 'Transmission_Coefficient'
    plot_parameters(title)

    # Find the file first!
    fileList = os.listdir(path)
    rexFileName = ''
    for file in fileList:
        if file[:18] == 'recrossing_factor_':
            rexFileName = file

    try:
        f = open(path + '/' + rexFileName, 'r')
    except FileNotFoundError:
        print('[ERROR] {} file not found! '.format(title))
    except PermissionError:
        print('[ERROR] {} file not found! '.format(title))
    else:
        fLines = f.readlines()
        f.close()
        time = []
        kappa = []
        for i in fLines[17:1500]:
            ele = i.split()
            time.append(np.float(ele[0]))
            kappa.append(np.float(ele[-1]))

        plt.xlabel('$t$ / fs')
        plt.ylabel('$\kappa(t)$')

        # plt.xscale('log')

        plt.xlim(time[0], time[-1])
        # endRF = np.mean(kappa[-5:])
        plt.axhline(y=kappa[-1], c=color[0], lw=0.5, linestyle='--')

        plt.plot(time, kappa, c=color[0])
        plt.tick_params(direction='in')
        # plt.show()
        plot_save(title)

if __name__=="__main__":
    path = input_path()
    xiList = get_xilist(path)
    plot_overlap(path, xiList)
    plot_variance(path, xiList)
    plot_pmf(path)
    plot_rexFactor(path)

