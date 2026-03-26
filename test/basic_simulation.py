import sys
import os
# import modules from root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))) 

from src.system import SolarSystem
import matplotlib.pyplot as plt

def main():
    """
    Runs simulation and animation using paraeterms of basic configuration
    """
    system = SolarSystem()
    system.read_parameters("base_parameters.json")
    system.run_simulation()
    _ = system.animation()

    plt.show()

if __name__ == "__main__":
    main()
    
