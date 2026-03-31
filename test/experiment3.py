import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# import modules from root
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, ROOT)
from src.system import SolarSystemSatelliteGrid

def heatmap(system):

    # plot heatmap of approaches
    fig, ax = system.plot_heatmap()
    ax.set_title("Initial Velocites' Closeness to Mars")
    plt.savefig(os.path.join(ROOT, 'figures', 'experiment3_heatmap.png'))
    plt.close()

def best_satellites(system):

    v_within = system.close_to_mars_satellites()
    print(f"We have {len(v_within)} satellites that get within 1e6 km of Mars (within 2 years):")
    for v in v_within:
        dist_mars, time_mars, _, _ = system.stats_from_velocity(v[0], v[1])
        print(f"v={v} gets within {dist_mars} AU of Mars after {time_mars:.4f} years")

    return v_within


def gets_close():
    # path to parameters
    file_name = "experiment3a.json"
    file_path = os.path.join(ROOT, 'parameters', file_name)

    system1 = SolarSystemSatelliteGrid()
    system1.read_parameters(file_path)

    vx_range = (0, 1, 100)
    vy_range = (6.3, 8.3, 100)

    system1.add_satellite_grid(vx_range, vy_range)
    system1.run_simulation()

    heatmap(system1)
    v_within = best_satellites(system1)

    return v_within

def return_to_earth(v_within):
    # path to parameters
    file_name = "experiment3b.json"
    file_path = os.path.join(ROOT, 'parameters', file_name)

    system2 = SolarSystemSatelliteGrid()
    system2.read_parameters(file_path)

    vx_range = (0, 1, 100)
    vy_range = (6.3, 8.3, 100)

    system2.add_satellite_grid(vx_range, vy_range)
    system2.run_simulation()

    print("and these satellites closest return to earth (within 10 years) are:")
    for v in v_within:
        _, _, dist_earth, time_earth = system2.stats_from_velocity(v[0], v[1])
        print(f"v={v} gets back to within {dist_earth} AU of Earth after {time_earth:.4f} years")


def main():
    v_within = gets_close()
    return_to_earth(v_within)


if __name__ == "__main__":
    main()

