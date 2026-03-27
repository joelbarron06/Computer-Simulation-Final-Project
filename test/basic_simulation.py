import sys
import os

# import modules from root
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

sys.path.insert(0, ROOT)
from src.system import SolarSystem

import matplotlib.pyplot as plt

def main():
    """
    Runs simulation and animation using paraeterms of basic configuration
    """
    system = SolarSystem()
    
    file_name = "base_parameters.json"
    file_path = os.path.join(ROOT, 'parameters', file_name)
    system.read_parameters(file_path)

    system.run_simulation()
    _ = system.animation()

    plt.show()

if __name__ == "__main__":
    main()
    
