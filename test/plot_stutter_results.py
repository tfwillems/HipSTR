import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

import math
import sys

def plot_dataset(datasets, title, pp):
    sort_keys   = [lambda x: (x[0][3:6], x[0][6:], x[0][0:3])]
    #sort_keys   = [lambda x: (x[0][4:7], x[0][0],  x[0][1], x[0][3], x[0][2])]

    orig_labels = [r"$log_{10}\ \mu$", r"$\beta$", r"$\rho_{m}$", r"$n_{samples}$" , r"$\rho_{s}$", r"$d$", r"$u$"]
    plot_vals   = [True, True, True, True, True, True, True]
    plot_dists  = [False, False, False, False, True, True, True]
    value_idx   = [-1,    -1,    -1,    -1, 0,    1,    2]

    for param in xrange(1):
        items      = sorted(datasets.items(), key = sort_keys[param])
        keys       = map(lambda x: x[0], items)
        fig, axes  = plt.subplots(sum(plot_vals), 1, sharex=True, sharey=False)
        positions  = range(len(datasets.keys()))
        subplot_id = 0
    
        for i in xrange(7):
            print(i)
            if not plot_vals[i]:
                continue

            if plot_dists[i]:
                val_sets = map(lambda x: map(lambda y: y[value_idx[i]], x[1]), items)
                boxes = axes[subplot_id].boxplot(val_sets, positions=positions, sym='')
                axes[subplot_id].set_ylabel(orig_labels[i])
            
                for line in boxes['whiskers']:
                    line.set_linestyle('-')

                prev_x_start = -0.25
                prev_key     = keys[0][i]
                for j in xrange(len(keys)):
                    if keys[j][i] != prev_key:
                        axes[subplot_id].plot([prev_x_start, (j-0.75)], [prev_key, prev_key], color='g')
                        prev_x_start = j-0.25
                        prev_key     = keys[j][i]
                axes[subplot_id].plot([prev_x_start, (len(keys)-0.75)], [prev_key, prev_key], color='g')
            else:
                vals = map(lambda x: x[0][i], items)
                axes[subplot_id].set_ylabel(orig_labels[i])

                pad = 0.25
                prev_x_start = -pad
                prev_val     = vals[0]
                for j in xrange(len(vals)):
                    if vals[j] != prev_val:
                        axes[subplot_id].plot([prev_x_start, (j-1+pad)], [prev_val, prev_val], color='g', linewidth=4, solid_capstyle='butt')
                        prev_x_start = j-pad
                        prev_val     = vals[j]
                axes[subplot_id].plot([prev_x_start, (len(keys)-1+pad)], [prev_val, prev_val], color='g', linewidth=4, solid_capstyle='butt')

            axes[subplot_id].tick_params(axis='both', which='major', labelsize=6)
            subplot_id += 1

        axes[0].set_title(title)
        axes[-1].xaxis.set_ticklabels([])
        axes[-1].set_xlabel("Simulation scenario")
        axes[0].set_ylim((-5.5,-1.0))
        axes[1].set_ylim((-0.05, 0.75))
        axes[2].set_ylim((0.55, 1.10))
        axes[3].set_ylim((50, 550))
        axes[4].set_ylim((0.5, 1.0))
        axes[5].set_ylim((-0.03, 0.5))
        axes[6].set_ylim((-0.04, 0.5))
        map(lambda x: x.xaxis.set_ticks_position('none'), axes)
        map(lambda x: x.yaxis.set_ticks_position('left'), axes)
        for ax in axes:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
        pp.savefig(fig)
        plt.close(fig)
    

def main():
    input_file  = sys.argv[1]
    output_file = sys.argv[2]
    data        = open(input_file, "r")
    datasets    = {}
    exc_count   = 0
    genotyper_indices = {}
    for line in data:
        if "EM_FAILED_TO_CONVERGE" in line:
            continue

        tokens = line.strip().split()

        mu, beta, p_geom        = map(float, tokens[0:3])
        pg_stutter,down,up      = map(float, tokens[3:6])
        est_pg,est_down,est_up  = map(float, tokens[10:13])

        mu = math.log10(mu)
        nsamps        = int(tokens[7])
        read_counts   = tokens[8]
        #read_percents = tokens[18]

        haploid    = tokens[6]
        phase_freq = tokens[9]

        key_one = (read_counts, haploid, phase_freq)
        if key_one not in datasets:
            datasets[key_one] = {}
        key_two = (mu, beta, p_geom, nsamps, pg_stutter, down, up)
        if key_two in datasets[key_one]:
            datasets[key_one][key_two].append([est_pg, est_down, est_up])
        else:
            datasets[key_one][key_two] = [[est_pg, est_down, est_up]]
    data.close()
    print("Excluding %d results to due to convergence failure messages"%(exc_count))

    pp = PdfPages(output_file)

    for key,dataset in sorted(datasets.items(), key=lambda x: (int(x[0][0].split(",")[0]), x[0][1:])):
        print(key)

        title = "Read Counts = " + key[0]
        if key[1] == "True":
            title += ", Haploid"
        else:
            title += ", Diploid"
            title += ", %% Phased Reads = %.1f"%(100.0*float(key[2]))


        if (key[1] == "True") and key[2] != "0.01":
            continue

        

        plot_dataset(dataset, title, pp)
    pp.close()    

if __name__ == "__main__":
    main()
