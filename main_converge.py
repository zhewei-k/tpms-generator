'''
------------
TPMS GENERATOR AND ANALYSIS TOOL v1.0.0
------------
Description: See README.md for guidance on how to operate.
Author: Zhe-Wei Kho
Last Updated: 12th June 2025
'''

import ntpms, digital_lab, logger, argparse, os
import matplotlib.pyplot as plt
import numpy as np

# @logger.timer #for visualisation
def main():

    #PARAMETERS
    global latticetype_0, cellsize_0, target_stiffness_0, mesh_available_0, plot_figure_0, custom_0
    latticetype_0 = 0       #global, constant
    cellsize_0 = 3.0        #global, constant mm
    thickness_0 = 0.1
    ext_force_0 = 1e3

    #SWITCHES
    mesh_available_0 = False
    custom_0 = r'C:\\Users\\User\\Documents\\GitHub\\fyp\\data\\2025_03_31_2107' ##if true
    custom_subfolder_0 = r'0'
    plot_figure_0 = True

    #CONVERGE EXPONENTIAL ONLY PARAMETERS
    start_iteration_0 = 3
    stop_iteration_0 = 5

    #CONVERGE CUSTOM ONLY PARAMETERS
    multiplier_list_0 = [
        # 2.0,
        1.0,
        # 0.5,
        # 1.0/3,
        # # 1.0/2,
        # 1.0/3
    ]

    #SECANT ONLY
    target_stiffness_0 = 3e9  #Pa
    max_iteration_0 = 5
    tolerance_0 = 0.05 #mm

    #argument parser
    parser = argparse.ArgumentParser(description="TPMS Analysis Argument Parser")
    parser.add_argument("mode", help="Analysis modes available: secant, converge, custom_converge") #required
    parser.add_argument("-l", "--lattice", type=int, help="0-Gyroid, 1-Schwarz, 2-Diamond, 3-Lidinoid, 4-SplitP, 5-Neovius") #optional
    args = parser.parse_args()
    if args.lattice:
        latticetype_0 = args.lattice
        print(f'### Lattice Argument Updated {latticetype_0}')
    print(f'### Lattice Argument Value is {latticetype_0}')

    print("### TPMS ANALYSIS STARTED")

    #folder management
    if mesh_available_0 == True:
        directory = logger.folder(custom_0)
        custom_0 = os.path.join(custom_0,custom_subfolder_0,'tpms_mesh.cdb' )
    elif mesh_available_0 == False:
        directory = logger.folder() #directory can be accessed via export.folder.directory

    #modes

    #CONVERGE EXPONENTIAL
    if args.mode == 'converge':
        multiplier_final_list, stiffness_final_list, ext_count_list = converge_exponential(latticetype_0, cellsize_0, thickness_0, start_iteration_0, stop_iteration_0, ext_force_0)
        print("### TPMS ANALYSIS COMPLETE")
        print(f"### Data has been stored in {directory}")
        for i in range(0,ext_count_list):
            print(f"### [{i}] --- MULTIPLIER (X) {multiplier_final_list[i]}")
            print(f"### [{i}] --- STIFFNESS (Pa) {stiffness_final_list[i]}")

    #CUSTOM CONVERGE
    elif args.mode == 'custom_converge':
        multiplier_final_list, stiffness_final_list, ext_count = converge_custom(latticetype_0, cellsize_0, thickness_0, multiplier_list_0, ext_force_0)
        print("### TPMS ANALYSIS COMPLETE")
        print(f"### Data has been stored in {directory}")
        for i in range(0,ext_count):
            print(f"### [{i}] --- MULTIPLIER (X) {multiplier_final_list[i]}")
            print(f"### [{i}] --- STIFFNESS (Pa) {stiffness_final_list[i]}")

    #SECANT
    elif args.mode == 'secant':
        thickness_final, stiffness_final, ext_count = secant(thickness_0, max_iteration_0, tolerance_0, ext_force_0)
        print("### TPMS ANALYSIS COMPLETE")
        print(f"### Data has been stored in {directory}")
        print(f"### [{ext_count}] --- THICKNESS (mm) {thickness_final}")
        print(f"### [{ext_count}] --- STIFFNESS (Pa) {stiffness_final}")
    
    
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
    def simulate_mesh_main(lat, cell, thick, elementmultiplier, force, count):
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
            distributed_force = np.divide(force, np.multiply(lengthnum, widthnum)) # divide force over number of unit cells, sigma = force / area
            print(f"### [{count}] --- DISTRIBUTED FORCE: {distributed_force}")
            stress, displacement = digital_lab.compressionV3(mesh, distributed_force, count, plot_figure_0)         # run simulation in Ansys
            strain = displacement / cell
            print(f'### [{count}] --- Z-AXIS STRAIN: {strain}')
            # heightnum = 5 #unit cells
            # distributed_strain = np.multiply(strain, heightnum) # divide strain over number of unit cells, strain = displacement / height
            
            return distributed_force, stress, strain
        
        mesh = "tpms_mesh.cdb"

        print(f"### [{count}] --- TPMS thickness {thick}")
        if mesh_available_0 == False:
            count, latticetype, cellsize, thickness, porosity, element_size = ntpms.mesh(lat, cell, thick, elementmultiplier, count)                   # generate mesh in nTop
            force, internal_stress, strain = extrapolate_cuboid(mesh, force, count)           #run simulation in Ansys as a test 
        elif mesh_available_0 == True:
            element_size, porosity = ntpms.mesh_stats(lat, cell, thick, elementmultiplier, count)
            force, internal_stress, strain = extrapolate_cuboid(custom_0, force, count)           #run simulation in Ansys as a test 

        # _, strain = digital_lab.compressionV2(mesh, force, count)         # run simulation in Ansys
        stress = force / (cellsize_0 * cellsize_0 * 1e-6)     # convert back to m^2
        stiffness = stress / strain

        logger.dataset['Count'] = count
        logger.dataset['LatticeType'] = latticetype 
        logger.dataset['CellSize'] = cellsize
        logger.dataset['Thickness'] = thickness
        logger.dataset['Porosity'] = porosity
        logger.dataset['ElementSize'] = element_size
        logger.dataset['Force'] = force
        logger.dataset['InternalStress'] = internal_stress
        logger.dataset['Strain'] = strain
        logger.dataset['Stiffness'] = stiffness

        return stiffness
        
    stiffness = simulate_mesh_main(lat, cell, thick, elementmultiplier, force, count)
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
        count (int)
    """    
    def define_secant(CSV_file, x, max_it, tol, force, count = None, x_prev = None, stiffness_prev = None):        
        if count is None:
            count = 0
            x_prev = np.add(x, np.multiply(tol,1.0005))   # to estimate x_prev
            stiffness_prev = simulate_mesh(CSV_file,latticetype_0, cellsize_0, x_prev, 1, force, count) # changed multiplier to 0.3
            # stiffness_prev = 646196877.7888234
        
        def char_equation(stiffness):
            return stiffness - target_stiffness_0

        try:
            count += 1
            stiffness = simulate_mesh(CSV_file,latticetype_0, cellsize_0, x, 1, force, count)

            grad = np.divide(
                np.subtract(char_equation(stiffness), char_equation(stiffness_prev)), #calculate new parameters
                np.subtract(x, x_prev)
            )

            x_new = np.subtract(x, np.divide(char_equation(stiffness), grad))
            print(f"### [{count}] --- grad: {grad}")
            print(f"### [{count}] --- x_new: {x_new}")

            if x_new <= 0:
                raise Exception("Next iteration thickness is a negative value.")

            if count >= max_it:         # analysis stops when we reach max_it
                print(f"### [{count}] ---  Exceeded maximum iterations.")
                print(f"### {count} Consider changing initial conditions (x, max_it or tol).")
                return x, stiffness, count

            if abs(char_equation(stiffness))/stiffness < 0.05:
                print(f"### [{count}] --- Optimisation complete.")
                print(f"### [{count}] --- x is {x}; x_new is {x_new}; x-x_new is {abs(x - x_new)}")
                print(f"### [{count}] --- x is {x}; x_new is {x_new}; x-x_new is {np.abs(np.subtract(x, x_new))}")
                return x, stiffness, count
            
            return define_secant(CSV_file, x_new, max_it, tol, force, count, x, stiffness)
        except Exception as e:
            print(f"### [{count}] --- An error occurred: {e}")
            quit()  #exit

    CSV_file = "tpms_secant.csv"
    # logger.spreadsheet(logger.dataset,CSV_file)

    x_final, stiffness_final, count = define_secant(CSV_file, x, max_it, tol, force)

    return x_final, stiffness_final, count
    
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
    nearest_stiffness = min(stiffness_list, key=lambda x: abs(x - target_stiffness_0))          # calculates absolute difference and find minimum value
    nearest_stiffness_id = stiffness_list.index(nearest_stiffness)                          # finds index of minimum value
    nearest_x = x_list[nearest_stiffness_id]    # finds corresponding thickness
    print(f"### [{count}] --- TARGET STIFFNESS: {target_stiffness_0/1e9} GPa")
    print(f"### [{count}] --- CLOSEST STIFFNESS: {nearest_x/1e9} GPa")

def converge_custom(lat, cell, thick, multiplier_list, force:float):

    count = 0
    CSV_file = "tpms_convergence_linear.csv"
    # logger.spreadsheet(logger.dataset,CSV_file)

    # large_multiplier = int(large_multiplier*1000)
    # small_multiplier = int(small_multiplier*1000)
    # step = int(step*1000)
    # multiplier_list = range(large_multiplier, small_multiplier, -step)
    # multiplier_list = [float(x) / 1000.0 for x in multiplier_list]
    print(multiplier_list)

    stiffness_list = []
    stress_list = []

    for i in multiplier_list:
        stiffness = simulate_mesh(CSV_file, lat, cell, thick, i, force, count)                 # generate mesh in nTop
        stress = logger.dataset['InternalStress']
        stress_list.append(stress/1e6)
        stiffness_list.append(stiffness)
        count +=1
    
    plt.style.use('_mpl-gallery')

    # data
    x = multiplier_list
    y = stress_list

    plt.figure(1, layout='constrained')
    plt.scatter(x, y, color='blue', label='Data points')
    plt.plot(x, y, color='red', label='Line')

    plt.xlabel('multiplier')
    plt.ylabel('Max Stress (MPa)')
    plt.title('Mesh Convergence with Exponential Increase in Mesh Size')
    plt.legend()
    plt.grid(True)

    return multiplier_list, stiffness_list, count

def converge_exponential(lat, cell, thick, start_val:float, stop_val:float, force:float):
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
    for i in range(start_val, stop_val):
        multiplier = 1.0 / float(2**i)
        multiplier_list.append(multiplier)

    stress_list = []
    stiffness_list = []
    for j in multiplier_list:
        stiffness = simulate_mesh(CSV_file, lat, cell, thick, j, force, count)                 # generate mesh in nTop
        stress = logger.dataset['InternalStress']

        stiffness_list.append(stiffness)
        stress_list.append(stress/1e6)        
        count +=1
    
    plt.style.use('_mpl-gallery')

    # data
    x = multiplier_list
    y = stress_list

    plt.figure(1, layout='constrained',figsize=(5, 3))
    plt.scatter(x, y, color='blue', label='')
    plt.plot(x, y, color='red', label='Line')

    plt.xlabel('multiplier')
    plt.ylabel('Stress (MPa)')
    plt.title('Mesh Convergence with Exponential Increase in Mesh Size')
    plt.grid(True)
    plt.savefig('mesh_convergence.png')

    plt.show()

    return multiplier_list, stiffness_list, count
    

if __name__ == "__main__":
    main()