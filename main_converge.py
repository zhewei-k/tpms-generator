'''
TPMS GENERATOR AND ANALYSIS TOOL
Author: Zhe-Wei Kho and Vee San Cheong
Description: Script that generates TPMS solids as STL files and CDB files and conducts simple analysis on them.
'''

import ntpms, digital_lab, logger
import matplotlib.pyplot as plt
import numpy as np

#initial values
latticetype_0 = 0       #global, constant
cellsize_0 = 3.0        #global, constant mm
target_stiffness = 5.5e9  #Pa

# @logger.timer #for visualisation
def main():
    print("### TPMS ANALYSIS STARTED")

    directory = logger.folder() #directory can be accessed via export.folder.directory
    thickness_final_list, stiffness_final_list, ext_count = converge_exponential(latticetype_0, cellsize_0, 0.5, 3, 1000.0)

    print("### TPMS ANALYSIS COMPLETE")
    print(f"Data has been stored in {directory}")
    for i in range(0,ext_count):
        print(f"[{i}] --- THICKNESS (mm) {thickness_final_list[i]}")
        print(f"[{i}] --- STIFFNESS (Pa) {stiffness_final_list[i]}")
    
    # os.chdir("..") #reset to original directory

def simulate_mesh(CSV_file, lat, cell, thick, elementmultiplier, force, count):
    """Generates a new mesh and runs a simulation for latest parameters

    Args:
        latticeType (str): 0-Gyroid, 1-Schwarz, 2-Diamond, 3-Lidinoid, 4-SplitP, 5-Neovius
        cellSize (float): size of unit cell (mm)
        thickness (float): thickness of tpms unit cell (mm)
        count (int) : folder reference number eg ./data/count

    Returns:
        summary (list): includes lattice type, cell size, thickness, element multiplier and count.
        stiffness (float): 
    """
    @logger.timer #this timer is logged!
    def extrapolate_cuboid(mesh, force, count):
        """Calculates forces for a 4x4x5 unit cell cuboid specimen

        Args:
            force (float): force applied on specimen (N)
            strain (float): displacement with respect to length of specimen in direction of force
            count (int): folder reference number eg ./data/count

        Returns:
            distributed_force (float): calculated force applied on each unit cell
            distributed_strain (float): calculated strain received by each unit cell
        """
        lengthnum = 1
        widthnum = 1
        heightnum = 5 #unit cells

        distributed_force = np.divide(force, np.multiply(lengthnum, widthnum)) # divide force over number of unit cells, sigma = force / area
        _, strain = digital_lab.compressionV3(mesh, distributed_force, count)         # run simulation in Ansys
        # distributed_strain = np.multiply(strain, heightnum) # divide strain over number of unit cells, strain = displacement / height
        
        return distributed_force, strain
    
    mesh = "tpms_mesh.cdb"

    print(f"### {count} TPMS thickness {thick}")
    count, latticetype, cellsize, thickness, porosity, element_size = ntpms.mesh(lat, cell, thick, elementmultiplier, count)                   # generate mesh in nTop
    # _, strain = digital_lab.compressionV2(mesh, force, count)         # run simulation in Ansys
    force, strain = extrapolate_cuboid(mesh, force, count)           #run simulation in Ansys as a test 
    stress = force / (cellsize_0 * cellsize_0 * 0.000001)     # convert back to m^2
    stiffness = stress / strain

    logger.dataset['Count'] = count
    logger.dataset['LatticeType'] = latticetype 
    logger.dataset['CellSize'] = cellsize
    logger.dataset['Thickness'] = thickness
    logger.dataset['Porosity'] = porosity
    logger.dataset['ElementSize'] = element_size
    logger.dataset['Force'] = force
    logger.dataset['Strain'] = strain
    logger.dataset['Stiffness'] = stiffness 
    logger.spreadsheet(logger.dataset,CSV_file) #writes to CSV

    print(logger.dataset.values())

    # summary.append(strain)
    # summary.append(force)
    # summary.append(stiffness)
    # output_file = open(csv_fname, mode = 'a')                           # append to csv
    # output_file.write(",".join([str(i) for i in summary]) + '\n')       # write to csv by inputting a piece of data and then creating a blank line
    # output_file.close()                                                 # close csv
    # print(summary)                                                      # [count, latticetype, cellsize, thickness, porosity, element_size, strain, force, stiffness]
    # print("\n")

    return stiffness

def secant(x, max_it, tol, force):
    """
    Secant method is a root-finding-method that optimises using thickness.

    Args:
        x (float): estimated value of the actual thickness
        max_it (int): _description_
        tol (float): to estimate the value of the previous thickness by a tolerance value
        force (float): force applied on specimen (N)
        count (int, optional): _description_. Defaults to None.
        x_prev (float, optional): _description_. Defaults to None.
        stiffness_prev (float, optional): _description_. Defaults to None.

    Raises:
        Exception: If next iteration thickness is a negative value

    Returns:
        x (float): optimised unit cell thickness (mm)
        stiffness (float): optimised unit cell stiffness matched with bone stiffness(Pa)
    """    
    def define_secant(CSV_file, x, max_it, tol, force, count = None, x_prev = None, stiffness_prev = None):        
        if count is None:
            count = 0
            x_prev = np.add(x, np.multiply(tol,1.0005))   # to estimate x_prev
            _, stiffness_prev = simulate_mesh(CSV_file,latticetype_0, cellsize_0, x_prev, 1, force, count) # changed multiplier to 0.3
            # stiffness_prev = 646196877.7888234
        
        def char_equation(stiffness):
            return stiffness - target_stiffness

        try:
            count += 1
            stiffness = simulate_mesh(CSV_file,latticetype_0, cellsize_0, x, 1, force, count)

            grad = np.divide(
                np.subtract(char_equation(stiffness), char_equation(stiffness_prev)), #calculate new parameters
                np.subtract(x, x_prev)
            )

            x_new = np.subtract(x, np.divide(char_equation(stiffness), grad))
            print(f"### {count} grad: {grad}")
            print(f"### {count} x_new: {x_new}")

            if x_new <= 0:
                raise Exception("Next iteration thickness is a negative value.")

            if count >= max_it:         # analysis stops when we reach max_it
                print(f"### {count} Exceeded maximum iterations.")
                print(f"### {count} Consider changing initial conditions (x, max_it or tol).")
                return x, stiffness

            if abs(char_equation(stiffness)) < 1e8:
                print(f"### {count} Optimisation complete.")
                print(f"x is {x}; x_new is {x_new}; x-x_new is {abs(x - x_new)}")
                print(f"x is {x}; x_new is {x_new}; x-x_new is {np.abs(np.subtract(x, x_new))}")
                return x, stiffness   
            
            return define_secant(CSV_file, x_new, max_it, tol, force, count, x, stiffness)
        except Exception as e:
            print(f"An error occurred: {e}")
            quit()  #exit

    CSV_file = "tpms_secant.csv"
    logger.spreadsheet(logger.dataset,CSV_file)

    x_final, stiffness_final = define_secant(CSV_file, x, max_it, tol, force)

    return x_final, stiffness_final
    
def process_list(max_x:float, min_x:float, step:float, force:float):
    """Processes a list of values from the maximum to the minimum for the steps specified

    Args:
        max_x (float): maximum value of range, inclusive (mm)
        min_x (float): minimum value of range, non-inclusive (mm)
        step (float): step size (mm)
        force (float): force applied on specimen (N)

    Returns:
        nearest_x: nearest corresponding thickness of nearest_stiffness (mm)
        nearest_stiffness: nearest stiffness of unit cell design to the stiffness of bone (Pa)
    """

    count = 0
    
    CSV_file = "tpms_list.csv"
    # logger.spreadsheet(logger.dataset,CSV_file)

    max_x = int(max_x*1000)
    min_x = int(min_x*1000)
    step = int(step*1000)
    x_list = range(max_x, min_x, -step)
    x_list = [float(x) / 1000.0 for x in x_list]
    print(x_list)

    stiffness_list = []

    for x in x_list:
        stiffness = simulate_mesh(CSV_file, latticetype_0, cellsize_0, x, 1, force, count)                 # generate mesh in nTop
        stiffness_list.append(stiffness)                                                    # run simulation in Ansys
        count +=1
    
    #finds closest correlation
    nearest_stiffness = min(stiffness_list, key=lambda x: abs(x - target_stiffness))          # calculates absolute difference and find minimum value
    nearest_stiffness_id = stiffness_list.index(nearest_stiffness)                          # finds index of minimum value
    nearest_x = x_list[nearest_stiffness_id]    # finds corresponding thickness
    print(f"Target stiffness: {target_stiffness/1e9} GPa")
    print(f"Closest stiffness: {nearest_x/1e9} GPa")

def converge_linear(lat, cell, thick, large_elementmultiplier:float, small_elementmultiplier:float, step:float, force:float):

    count = 0
    
    CSV_file = "tpms_convergence_linear.csv"
    # logger.spreadsheet(logger.dataset,CSV_file)

    large_elementmultiplier = int(large_elementmultiplier*1000)
    small_elementmultiplier = int(small_elementmultiplier*1000)
    step = int(step*1000)
    multiplier_list = range(large_elementmultiplier, small_elementmultiplier, -step)
    multiplier_list = [float(x) / 1000.0 for x in multiplier_list]
    print(multiplier_list)

    stiffness_list = []

    for i in multiplier_list:
        stiffness = simulate_mesh(CSV_file, lat, cell, thick, i, force, count)                 # generate mesh in nTop
        stiffness_list.append(stiffness)
        count +=1
    
    plt.style.use('_mpl-gallery')

    # data
    x = multiplier_list
    y = stiffness_list

    plt.figure(1, layout='constrained')
    plt.scatter(x, y, color='blue', label='Data points')
    plt.plot(x, y, color='red', label='Line')

    plt.xlabel('multiplier')
    plt.ylabel('stiffness (Pa)')
    plt.title('Mesh Convergence with Exponential Increase in Mesh Size')
    plt.legend()
    plt.grid(True)

    return multiplier_list, stiffness_list, count

def converge_exponential(lat, cell, thick, stop_val:float, force:float):
    """Conducts a mesh convergence check to investigate the mesh convergence capabilities

    Args:
        latticeType (str): 0-Gyroid, 1-Schwarz, 2-Diamond, 3-Lidinoid, 4-SplitP, 5-Neovius
        cellSize (float): size of unit cell (mm)
        thickness (float): thickness of tpms unit cell (mm)
        count (int) : folder reference number eg ./data/count

    Returns:
        multiplier_list, stiffness_list, count
    """
    count = 0
    
    CSV_file = "tpms_converge_exponential.csv"
    # logger.spreadsheet(logger.dataset,CSV_file)

    multiplier_list = []
    for i in range(0, stop_val):
        multiplier = 1.0 / float(2**i)
        multiplier_list.append(multiplier)

    stiffness_list = []
    for j in multiplier_list:
        stiffness = simulate_mesh(CSV_file, lat, cell, thick, j, force, count)                 # generate mesh in nTop
        stiffness_list.append(stiffness)
        count +=1
    
    plt.style.use('_mpl-gallery')

    # data
    x = multiplier_list
    y = stiffness_list

    plt.figure(1, layout='constrained')
    plt.scatter(x, y, color='blue', label='Data points')
    plt.plot(x, y, color='red', label='Line')

    plt.xlabel('multiplier')
    plt.ylabel('stiffness (Pa)')
    plt.title('Mesh Convergence with Exponential Increase in Mesh Size')
    plt.legend()
    plt.grid(True)

    plt.show()

    return multiplier_list, stiffness_list, count
    

if __name__ == "__main__":
    main()