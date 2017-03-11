#-*- coding: utf-8 -*-
"""
Created on Wed Feb 28 2017
@author: patrick.diehl@polymtl.ca
@author: ilyass.tabiai@polymtl.ca
@author: rolland.delorme@polymtl.ca
"""



import matplotlib.pyplot as plt

## Class handles plotting of results using Matplorlib, exports plots as PDFs
class Plotting():
     ## Plots energy at a certain time step for each nodes, currentyl only for 1D results
     # @param energy Array containing energy values for each node at each time step
     # @param time Moment at which results should be plotted
     # @param initial 
     # @param outpath Desired path for the output PDF
     
     def plot_energy(self, energy, time, initial, outpath):
        # print len(energy) , len(time), len(initial)
        maxvalues = []
        color = []
        for i in range(0, len(initial)):
            e = []
            for j in range(0, len(energy)):
                e.append(abs(energy[j][i]))
            line = plt.plot(time, e, marker="o")
            color.append(line[0].get_color())
            maxvalues.append(max(e))
        # print color
        # plt.plot([0, max(time)], [max(maxvalues) + 100,
        #                         max(maxvalues) + 100], lw=2, c='black')
        # for i in range(len(initial)):
        #    plt.plot(
        #        (max(time) / 37.5) * initial[i],
        #        max(maxvalues) + 100,
        #        color[i],
        #        marker='o')
        plt.grid()
        plt.xlabel("Time [s]")
        plt.ylabel("Strain Energy")
        plt.savefig(outpath +
                    "energy_" +
                    str(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")) +
                    ".pdf")

    #def plot_force(self, PD_deck):
    #    for t_n in range(1, PD_deck.Num_TimeStep):
    #        force_plot = plt
    #        force_plot.plot(self.y[:, 1], self.forces[:, t_n], '-+', label=t_n)
    #    force_plot.legend(title="forces")
    #    return force_plot

    #def plot_positions(self, PD_deck):
    #    for t_n in range(1, PD_deck.Num_TimeStep):
    #        position_plot = plt
    #        position_plot.plot(self.y[:, 1], self.y[:, t_n], '-+')
    #    position_plot.legend(title="position")
    #    return position_plot

