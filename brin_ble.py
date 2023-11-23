import glfw
import numpy as np

from scenes import *
from viewer import Viewer

	

        


# Separates the main to make sure all objects are deleted
# before glfw.terminate is called
def main():
    viewer=Viewer(height=960, width=1280,
                  bgColor = np.array([0.4, 0.4, 0.4]))

    # Loading the scene
    SystemPenduleTest(viewer, "implicit")
    
    # Main loop
    viewer.run()



if __name__ == '__main__':
	
    import sys
    
    # Initialization
    glfw.init()
	
	
    main()
    
    # End
    glfw.terminate()


