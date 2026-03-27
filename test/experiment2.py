import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

# import modules from root
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from src.system import SolarSystem, SolarSystemDirectEuler, SolarSystemEulerCromer

def simulate_data():
    """
    Simulate solar system using each integration algorithm.
    Write each energy to a unique file.
    """
    INTERVAL = 20 # with dt=0.01 gives 5 readings per year

    # path to parameters
    file_name = f"experiment2.json"
    file_path = os.path.join(ROOT, 'parameters', file_name)

    # beeman integration
    system_bmn = SolarSystem()
    system_bmn.read_parameters(file_path)
    system_bmn.run_simulation(print_periods=False, energy_file_name='e2_beeman_energy', energy_interval=INTERVAL)
    

    # euler cromer integration
    system_elcr = SolarSystemEulerCromer()
    system_elcr.read_parameters(file_path)
    system_elcr.run_simulation(print_periods=False, energy_file_name='e2_euler_cromer_energy', energy_interval=INTERVAL)
    

    # direct euler integration
    system_eul = SolarSystemDirectEuler()
    system_eul.read_parameters(file_path)
    system_eul.run_simulation(print_periods=False, energy_file_name='e2_direct_euler_energy', energy_interval=INTERVAL)
    
    

def combine_data():
    """
    Combines energy data of each method into a single time series data frame.
    """
    # read simulated csv files
    bmn_path = os.path.join(ROOT, 'data', 'e2_beeman_energy.csv')
    df_bmn = pd.read_csv(bmn_path)

    elcr_path = os.path.join(ROOT, 'data', 'e2_euler_cromer_energy.csv')
    df_elcr = pd.read_csv(elcr_path)

    eul_path = os.path.join(ROOT, 'data', 'e2_direct_euler_energy.csv')
    df_eul = pd.read_csv(eul_path)

    # join into one timeseries data frame
    df = pd.concat([df_bmn, df_elcr["energy"], df_eul["energy"]], axis=1)
    df.columns = ['time', 'Beeman', 'Euler-Cromer', 'Direct Euler']
    df = df.set_index("time")

    return df

def rolling_mean_plot(data):
    """
    Create timeseries plots of energy in each method with a rolling mean.
    """
    fig, ax = plt.subplots(3, figsize=(12, 10))

    plot_data = data.iloc[5:] # skip first year

    ax[0].plot(plot_data["Beeman"], label='System Energy')
    ax[0].plot(plot_data["Beeman"].rolling(window=50).mean(), 'r-', alpha=0.6, label='10yr Rolling Mean')
    ax[0].set(title='Beeman', xlabel='time (yr)', ylabel='Energy')
    ax[0].legend()

    ax[1].plot(plot_data["Euler-Cromer"], label='System Energy')
    ax[1].plot(plot_data["Euler-Cromer"].rolling(window=50).mean(), 'r-', alpha=0.6, label='10yr Rolling Mean')
    ax[1].set(title='Euler-Cromer', xlabel='time (yr)', ylabel='Energy')
    ax[1].legend()

    ax[2].plot(plot_data["Direct Euler"], label='System Energy')
    ax[2].plot(plot_data["Direct Euler"].rolling(window=50).mean(), 'r-', alpha=0.6, label='10yr Rolling Mean')
    ax[2].set(title='Direct Euler', xlabel='time (yr)', ylabel='Energy')
    ax[2].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, 'figures', 'rolling_mean_comparison'))
    plt.close()


def normalised_plot(data):
    """
    Plot normalised energy for each method on same axes.
    """
    plot_data = data.iloc[5:]
    
    for col in plot_data.columns:
        normalised = (plot_data[col] - plot_data[col].iloc[0])/abs(plot_data[col].iloc[0])
        plot_data[f"{col} Norm"] = normalised

    fig, ax = plt.subplots(2, figsize=(8,8))

    ax[0].plot(plot_data['Beeman Norm'], label='Beeman')
    ax[0].plot(plot_data['Euler-Cromer Norm'], label='Euler-Cromer')
    ax[0].plot(plot_data['Direct Euler Norm'], label='Direct Euler')
    ax[0].set(title='Normalised Energy for each Integration Method',  ylabel='Normalised Energy')
    ax[0].legend()

    ax[1].plot(plot_data['Beeman Norm'], label='Beeman', alpha=0.6)
    ax[1].plot(plot_data['Euler-Cromer Norm'], label='Euler-Cromer', alpha=0.6)
    ax[1].set(xlabel='time (yr)', ylabel='Normalised Energy')
    ax[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, 'figures', 'normalised_energy_comparison'))
    plt.close()

def stats_for_discussion(data):
    """
    Create csv file with statistics for use in discussion.
    """
    df = data.iloc[5:]

    with open(os.path.join(ROOT, "data", "e3_statistics.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['algorithm', 'mean', 'std', 'range', 'slope', 'monotone increasing']) # header

        for col_name in df.columns:
            mean = df[col_name].mean()
            std = df[col_name].std()
            monotone_inc = df[col_name].is_monotonic_increasing
            data_range = abs(df[col_name].min() - df[col_name].max())

            # linear regression; show drift in energy
            m, _ = np.polyfit(df.index, df[col_name], 1)

            writer.writerow([col_name, mean, std, data_range, m, monotone_inc])



def main():
    simulate_data()
    data = combine_data()
    rolling_mean_plot(data)
    normalised_plot(data)
    stats_for_discussion(data)


if __name__ == "__main__":
    main()