"""
------------
DIGITAL_LAB v1.0.0
------------
DIGITAL_LAB is a module that provides you with tools to perform analysis using code from custom meshes.
Requires the following files to be in the same folder directory
- input_template.json
- tpms.ntop
"""

import numpy as np
#import pyvista as pv
import os, psutil, logger
import ansys.mapdl.core as pymapdl

def main():
    custom = r'C:\\Users\\User\\Documents\\GitHub\\fyp\\data\\2025_05_06_1814'
    logger.folder(custom)
    custom = os.path.join(custom,'99','tpms_mesh.cdb' )
    compressionV3(custom, 100, 0, True) #remember to change the file location here too
    # main_fetch_enum()

def main_fetch_enum():
    custom = r'C:\\Users\\User\\Documents\\GitHub\\fyp\\data-final\\2025_03_31_2036 - 0.4'
    logger.folder(custom)
    count = '3'
    custom = os.path.join(custom,count,'tpms_mesh.cdb')
    fetch_enum(custom,count)

def find_process_pid(process_name):
    pid = None
    for proc in psutil.process_iter(['pid', 'name']):
        if proc.info['name'] == process_name:
            pid = proc.info['pid']
            break
    return pid

def fetch_enum(mesh,count):
    try:
        directory = logger.folder.directory
        sub_directory = os.path.join(directory,str(count))
        status = True
        mapdl = pymapdl.launch_mapdl(
            #jobname="Database Test",  #causes errors if turned on
            nproc=4,
            # loglevel="INFO",
            start_timeout=120.0,
            cleanup_on_exit=True,
            # log_apdl=os.path.join(sub_directory,r"pymapdl_compressionV3.txt"),
            additional_switches='-smp',
            # port=5000,
            exec_file=r"C:\\Program Files\ANSYS Inc\\v242\\ansys\\bin\\winx64\\ANSYS242.exe", # version 24.2
            start_instance=status,
            remove_temp_dir_on_exit=True
            # version=242
        )
        mapdl.clear()

        mesh_path = os.path.join(sub_directory,mesh)
        print(f"### [{count}] --- MESH PATH: "+mesh_path)
        mapdl.input(mesh_path)
        print(mapdl.cmlist())

        mapdl.prep7()
        mapdl.units("SI", mute=False)

        element_num = np.max(mapdl.esel('ALL'))
        logger.dataset['ElementNum'] = element_num
        print(f'\n### [{count}] --- ELEMENT NUMBER: {element_num}')

        mapdl.finish()
        mapdl.wait(1)
        mapdl.exit()

    except Exception as e:
        print(f"### [{count}] --- MAPDL LAUNCH ERROR ENCOUNTERED: {e}")
        try:
            mapdl.finish()
            mapdl.wait(1)
            mapdl.exit()
        except:
            print(f'### [{count}] --- MAPDL EXIT ERROR ENCOUNTERED')

        # I have turned off the termination function to avoid interfering with any other software that is running
        # will quit() after encountering an error
        quit()

    
@logger.timer #for visualisation
def compressionV3(mesh,force,count,plot_figures=True):

    directory = logger.folder.directory
    sub_directory = os.path.join(directory,str(count))
    # managing errors from launching PyMAPDL
    try:
        if count == 0:
            status = True
        else:
            status = False
        mapdl = pymapdl.launch_mapdl(
            #jobname="Database Test",  #causes errors if turned on
            nproc=4,
            # loglevel="INFO",
            start_timeout=120.0,
            cleanup_on_exit=True,
            # log_apdl=os.path.join(sub_directory,r"pymapdl_compressionV3.txt"),
            additional_switches='-smp',
            # port=5000,
            exec_file=r"C:\Program Files\ANSYS Inc\v242\ansys\bin\winx64\ANSYS242.exe", # version 24.2
            start_instance=status,
            remove_temp_dir_on_exit=True
            # version=242
        )
        print(f"### [{count}] --- MAPDL LAUNCH SUCCESS")
        print(f"### [{count}] --- SIMULATION STARTING")

        mapdl.clear()

        mesh_path = os.path.join(sub_directory,mesh)
        print(f"### [{count}] --- MESH PATH: "+mesh_path)
        mapdl.input(mesh_path)
        print(mapdl.cmlist())

        mapdl.prep7()
        mapdl.units("SI", mute=False)
        print(mapdl.etlist())

        # Redefine element type - hidden
        # element_type = "SOLID186"
        # mapdl.et(2,element_type)
        # mapdl.emodif("ALL", "TYPE", 2)
        # print(f"### [{count}] --- ELEMENT TYPE CHANGED TO {element_type}")

        # Check element type - hidden
        # elem_id = 1
        # elem_type = mapdl.get_value("ELEM", elem_id, "ATTR", "TYPE")
        # print(f"Element {elem_id} uses element {elem_type}")
        # mapdl.wait(5)

        # Define and assign material
        mapdl.mp("EX", 2, 110e9)  # Elastic moduli in Pa (kg/(m*s**2))
        mapdl.mp("DENS", 2, 4510)  # Density in kg/m3
        mapdl.mp("NUXY", 2, 0.34)  # Poisson's Ratio
        mapdl.emodif('ALL', 'MAT', 2)

        # Print element plot
        print(f"### [{count}] --- EPLOT STARTING")
        eplot_path = os.path.join(sub_directory,"eplot.jpg")
        mapdl.eplot(
            vtk=True, 
            show_edges=True,
            line_width=2, 
            background="w",
            savefig=eplot_path
        )
        print(f"### [{count}] --- EPLOT COMPLETE")

        # Apply constraint boundary condition
        base_nodes = mapdl.components["NODESET2"].items
        print(f"### [{count}] --- BASE NODES: {len(base_nodes)}")
        mapdl.d("NODESET2", "ALL", 0)

        # Apply force boundary condition (optionally, you can turn on displacement)
        # mapdl.d("NODESET1","UX",0)
        # mapdl.d("NODESET1","UY",0)
        # mapdl.d("NODESET1","UZ",-0.00001)
        mapdl.cp(1, "UZ", "NODESET1") #reference num, degrees of freedom, which nodes
        force_nodes = mapdl.components["NODESET1"].items
        print(f"### [{count}] --- FORCE NODES: {len(force_nodes)}")
        nodal_force = -force / len(force_nodes)
        mapdl.f("NODESET1", "FZ", nodal_force)

        # Calculate element number and save
        element_num = np.max(mapdl.esel('ALL'))
        node_num = np.max(mapdl.nsel('ALL'))
        logger.dataset['ElementNum'] = element_num
        logger.dataset['NodeNum'] = node_num
        print(f'\n### [{count}] --- ELEMENT NUMBER: {element_num}')
        print(f'\n### [{count}] --- NODE NUMBER: {node_num}')

        # Solve
        mapdl.allsel(mute=True)
        mapdl.run("/SOLU")
        mapdl.antype("STATIC")
        mapdl.wait(1)

        print(f'\n### [{count}] --- SOLVER STARTING')
        
        attempts = 0
        max_attempts = 1
        success = False
        while attempts < max_attempts or success != True:
            try:
                print(f'### [{count}] --- SOLVER ATTEMPT {attempts+1}')
                mapdl.solve(mute=True)
                success = True
            except Exception as e:
                print(f'\n### [{count}] --- SOLVER FAILED: {e}')
                print(f'\n### [{count}] --- CHECKING WHETHER MAPDL IS ALIVE')
                print(f'### [{count}] --- MAPDL STATUS: {mapdl.is_alive()}')
            attempts +=1
                
        mapdl.finish()
        print(f'\n### [{count}] --- SOLVER COMPLETE')

        mapdl.post1()
        mapdl.set(1)

        max_strain = np.nanmax(mapdl.post_processing.nodal_total_eqv_strain())
        print(f'### [{count}] --- MAX STRAIN-EQV: {max_strain}')
        
        max_stress = np.nanmax(mapdl.post_processing.nodal_eqv_stress())
        print(f'### [{count}] --- MAX STRESS: {max_stress/1e6} MPa')

        max_z_displacement = np.nanmax(np.abs(mapdl.post_processing.nodal_displacement('Z'))) #Z-direction displacement in meters
        max_z_displacement = max_z_displacement * 1e3
        print(f'### [{count}] --- MAX DISPLACEMENT: {max_z_displacement} mm')

        # alternative plotting method
        # mapdl.result.plot_nodal_displacement(
        #     rnum=0,
        #     show_displacement = True,
        #     displacement_factor = 20.0,
        # )

        if plot_figures == True:
            strain_plot_path = os.path.join(sub_directory,"strain_plot.jpg")
            mapdl.post_processing.plot_nodal_total_eqv_strain(
                cpos="iso",
                lighting=False,
                background="white",
                edge_color="black",
                show_edges=True,
                # n_colors=9,
                cmap="jet",
                font_size="20",
                # return_plotter=True,
                savefig=strain_plot_path,
            )
            stress_plot_path = os.path.join(sub_directory,"stress_plot.jpg")
            mapdl.post_processing.plot_nodal_eqv_stress(
                cpos="iso",
                lighting=False,
                background="white",
                edge_color="black",
                show_edges=True,
                # n_colors=9,
                cmap="jet",
                # return_plotter=True,
                savefig=stress_plot_path,
            )
            displacement_plot_path = os.path.join(sub_directory,"displacement_plot.jpg")
            mapdl.post_processing.plot_nodal_displacement(
                'Z',
                cpos="iso",
                lighting=False,
                background="white",
                edge_color="black",
                show_edges=True,
                # n_colors=9,
                cmap="jet",
                # cmap_min="-3e-7",
                # return_plotter=True,
                savefig=displacement_plot_path,
            )

        mapdl.finish()
        print(f"### [{count}] --- SIMULATION COMPLETE")
        return max_stress, max_z_displacement
    except Exception as e:
        # quit after encountering an error
        print(f"### [{count}] --- MAPDL LAUNCH ERROR ENCOUNTERED: {e}")
        try:
            mapdl.finish()
            mapdl.wait(1)
            mapdl.exit()
        except:
            print(f'### [{count}] --- MAPDL EXIT ERROR ENCOUNTERED')
        quit()  # quit after encountering an error

    
if __name__ == "__main__":
    main()
