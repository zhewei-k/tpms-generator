import ntpms, digital_lab, os, shutil
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

#initial values
latticetype_0 = 0
cellsize_0 = 3.0                      #mm
stiffness_bone = 5.5e9               #Pa

now = datetime.now()

def char_equation(stiffness):
    return stiffness - stiffness_bone

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
    lengthnum = 3
    widthnum = 3
    heightnum = 5 #unit cells

    distributed_force = np.divide(force, np.multiply(lengthnum, widthnum)) # divide force over number of unit cells, sigma = force / area
    _, strain = digital_lab.compressionV2(mesh, distributed_force, count)         # run simulation in Ansys
    # distributed_strain = np.multiply(strain, heightnum) # divide strain over number of unit cells, strain = displacement / height
    
    return distributed_force, strain

def simulate_mesh(csv_fname, lat, cell, thick, elementmultiplier, force, count):
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
    
    mesh = "tpms_mesh.cdb"

    print(f"### {count} TPMS thickness {thick}")
    summary = ntpms.mesh(lat, cell, thick, elementmultiplier, count)                   # generate mesh in nTop
    # _, strain = digital_lab.compressionV2(mesh, force, count)         # run simulation in Ansys
    force, strain = extrapolate_cuboid(mesh, force, count)           #run simulation in Ansys as a test 
    stress = force / (cellsize_0 * cellsize_0 * 0.000001)     # convert back to m^2
    stiffness = stress / strain

    summary.append(strain)
    summary.append(force)
    summary.append(stiffness)
    output_file = open(csv_fname, mode = 'a')                           # append to csv
    output_file.write(",".join([str(i) for i in summary]) + '\n')       # write to csv by inputting a piece of data and then creating a blank line
    output_file.close()                                                 # close csv
    print(summary)                                                      # [count, latticetype, cellsize, thickness, porosity, element_size, strain, force, stiffness]
    print("\n")

    return summary, stiffness

def define_secant(csv_fname, x, max_it, tol, force, count = None, x_prev = None, stiffness_prev = None):
    """Secant method is a root-finding-method that optimises using thickness.

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
    
    if count is None:
        count = 0
        x_prev = np.add(x, np.multiply(tol,1.0005))   # to estimate x_prev
        _, stiffness_prev = simulate_mesh(csv_fname,latticetype_0, cellsize_0, x_prev, 1, force, count) # changed multiplier to 0.3
        # stiffness_prev = 646196877.7888234
        
    try:
        count += 1
        _, stiffness = simulate_mesh(csv_fname,latticetype_0, cellsize_0, x, 1, force, count)

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
        
        return define_secant(csv_fname, x_new, max_it, tol, force, count, x, stiffness)
    except Exception as e:
        print(f"An error occurred: {e}")
        quit()  #exit

def execute_secant(x, max_it, tol, force):
    CSV_file = "tpms_secant_"+now.strftime("%Y_%m_%d_%H%M")+".csv"
    Header = ['Count', 'Lattice Type', 'Cell Size', 'Thickness', 'Porosity', 'Mesh Size', 'Strain', 'Force', 'Stiffness']
    output_file = open(CSV_file, mode='w')
    output_file.write(",".join(Header)+'\n')
    output_file.close()

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
    
    CSV_file = "tpms_list_"+now.strftime("%Y_%m_%d_%H%M")+".csv"
    Header = ['Count', 'Lattice Type', 'Cell Size', 'Thickness', 'Porosity', 'Mesh Size', 'Stress', 'Strain', 'Stiffness']
    output_file = open(CSV_file, mode='w')
    output_file.write(",".join(Header)+'\n')
    output_file.close()

    max_x = int(max_x*1000)
    min_x = int(min_x*1000)
    step = int(step*1000)
    x_list = range(max_x, min_x, -step)
    x_list = [float(x) / 1000.0 for x in x_list]
    print(x_list)

    stiffness_list = []

    for x in x_list:
        _, stiffness = simulate_mesh(CSV_file, latticetype_0, cellsize_0, x, 1, force, count)                 # generate mesh in nTop
        stiffness_list.append(stiffness)                                                    # run simulation in Ansys
        count +=1
    
    nearest_stiffness = min(stiffness_list, key=lambda x: abs(x - stiffness_bone))          # calculates absolute difference and find minimum value
    nearest_stiffness_id = stiffness_list.index(nearest_stiffness)                          # finds index of minimum value
    nearest_x = x_list[nearest_stiffness_id]    # finds corresponding thickness

    return nearest_x, nearest_stiffness

def converge_stiffness(lat, cell, thick, large_elementmultiplier:float, small_elementmultiplier:float, step:float, force:float):
    count = 0
    
    CSV_file = "tpms_convergence_"+now.strftime("%Y_%m_%d_%H%M")+".csv"
    Header = ['Count', 'Lattice Type', 'Cell Size', 'Thickness', 'Porosity', 'Mesh Size', 'Strain', 'Force', 'Stiffness']
    output_file = open(CSV_file, mode='w')
    output_file.write(",".join(Header)+'\n')
    output_file.close()

    large_elementmultiplier = int(large_elementmultiplier*1000)
    small_elementmultiplier = int(small_elementmultiplier*1000)
    step = int(step*1000)
    multiplier_list = range(large_elementmultiplier, small_elementmultiplier, -step)
    multiplier_list = [float(x) / 1000.0 for x in multiplier_list]
    print(multiplier_list)

    stiffness_list = []

    for i in multiplier_list:
        _, stiffness = simulate_mesh(CSV_file, lat, cell, thick, i, force, count)                 # generate mesh in nTop
        stiffness_list.append(stiffness)
                                                            # run simulation in Ansys
        count +=1
    
    plt.style.use('_mpl-gallery')

    # data
    x = multiplier_list
    y = stiffness_list

    plt.figure(1, layout='constrained')
    plt.scatter(x,y)
    plt.xlabel('multiplier')
    plt.ylabel('stiffness (Pa)')
    plt.show()

    return multiplier_list, stiffness_list, count

def main():
    print("---- TPMS ANALYSIS STARTED ----")

    if os.path.isdir('./data/'):    #checks whether "data" folder is present,
        shutil.rmtree('./data/')    #if yes, deletes "data" folder
    os.mkdir('./data/')             #creates folder if data folder for storing data
    
    thickness_final_list, stiffness_final_list, ext_count = converge_stiffness(0, 3, 0.66, 0.9, 0.2, 0.1, 100.0)

    print("---- TPMS ANALYSIS COMPLETE ----")
    print("Data has been stored in ./data/")
    for i in range(0,ext_count):
        print(f"THICKNESS (mm) {thickness_final_list[i]} ---  STIFFNESS (Pa) {stiffness_final_list[i]}")
    

if __name__ == "__main__":
    main()