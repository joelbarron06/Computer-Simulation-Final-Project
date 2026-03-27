import sys
import os

# import modules from root
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

sys.path.insert(0, ROOT)
from src.system import SolarSystem

def orbtial_periods(dt):
        system = SolarSystem()

        file_name = f"experiment1.json"
        file_path = os.path.join(ROOT, 'parameters', file_name)
        system.read_parameters(file_path)

        # overwrite timestep for function passed
        system.timestep = dt

        print("="*10, f"Orbital periods for dt = {system.timestep}:", sep="\n")

        system.run_simulation(write_energy=False)


def main():
      dt_vals = [0.1, 0.01, 0.001, 0.0005]

      for dt in dt_vals:
            orbtial_periods(dt)

if __name__ == "__main__":
    main()

