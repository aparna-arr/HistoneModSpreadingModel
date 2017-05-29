import Constants
import random
import numpy as np
import operator
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.stats as stats
import math
from MyEnum import ProbSpread, States, Domain, DomainBleed, Divisions, ProbConv

def animate_from_file(filename, n, f, t, outfile):
    global colors
    global totals
    colors = [Constants.GRAY]*n
    totals = {
        States.A_STATE : 0,
        States.M_STATE : 0,
        States.U_STATE : 0
            }
    # base figure
    fig = plt.figure()
    # axes
    # figure 1: nucleosomes
    ax1 = fig.add_subplot(2,2,1)
    ax1.set_xticks([])
    ax1.set_yticks([])

    titlestr = "N = " + str(n) + " F = " + str(f)
    ax1.set_title(titlestr)

    # figure 2: proportion vs events
    ax2 = fig.add_subplot(2,2,2)
    ax2.set_ylim([-5,105])
    ax2.set_xlim(t / 40, t + t / 40)
    ax2.set_xlabel("Timesteps")
    ax2.set_ylabel("% Nucleosomes")

#    ax3 = fig.add_subplot(2,2,3)
    # adjust spacing between figures
    fig.tight_layout()

    # calculate size of square for nucleosomes in fig 1
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    x_int = 1 / cols
    y_int = 1 / rows

    x_vals = []
    y_vals = []

    count = 0
    for i in range(rows):
        for j in range(cols):
            if count >= n :
                break
            
            x_vals.append(j * x_int)
            y_vals.append(i * y_int)
            count += 1

    # initialize fig1 as a scatterplot
    scat = ax1.scatter(x_vals, y_vals, facecolors = colors)

    # initialize fig2 as a lineplot
    lines = []
    lines.append(ax2.plot([], [], lw = 2, color = "red", label = "A")[0])
    lines.append(ax2.plot([], [], lw = 2, color = "blue", label = "M")[0])
    ax2.legend(loc = "upper right")

    for line in lines:
        line.set_data([],[])

    linex = []
    liney = [[], []]

    global file_sim
    file_sim = open(filename, "r")

    def init_an():
        '''initialize animation: this is needed'''
        return scat, (*lines)

    def update_an(i):
        '''update function for animation'''
        if i == 0:
            return scat, (*lines)

        line = file_sim.readline()

       # print("i is ",i)

        totals[States.A_STATE] = 0
        totals[States.U_STATE] = 0
        totals[States.M_STATE] = 0
        index = 0

        for nuc in list(line.rstrip()):
            #print("len colors:", len(colors), "index:", index)
            totals[States.string_to_enum(nuc)] += 1
            colors[index] = Constants.state_to_color(States.string_to_enum(nuc))
            index += 1

#            print("PASSING self.TIME as t:", self.TIME)
#            self.TIME = self.run_event(self.TIME)
#            self.time_array.append(self.TIME)
        #print("Frame:",i)
        #print("Time is:", self.TIME)
#            print("After event, self.TIME is", self.TIME)

#            if self.dat['d'] != Divisions.NONE and i != 0:
#                div = self.dat['data']['divisions']
#                if i % int(self.dat['e'] / div) == 0:
#                    self.divide()
#                    lines.append(ax2.plot([], [], lw = 1, ls = "dotted", color = "black")[0])
#                    lines[len(lines) - 1].set_data([i, i], [-5, 105])

        linex.append(i)
        liney[0].append(totals[States.A_STATE] / n * 100)
        liney[1].append(totals[States.M_STATE] / n * 100)

        lines[0].set_data(linex, liney[0])
        lines[1].set_data(linex, liney[1])

        scat.set_facecolors(colors)

        return scat, (*lines)

    anim = animation.FuncAnimation(fig, update_an, init_func = init_an, frames = t, interval = 1, repeat = False, blit=True)

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps = 200, metadata = dict(artist = 'Me'), bitrate = 1800)
    anim.save(outfile + ".mp4", writer = writer)
    file_sim.close()

