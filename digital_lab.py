"""
------------
DIGITAL_LAB v1.0.0
------------
DIGITAL_LAB is a module that provides you with tools to perform analysis using code from custom meshes.
Requires the following files to be in the same folder directory
- input_template.json
- tpms.ntop
"""

#from matplotlib import pyplot as plt
import numpy as np
#import pyvista as pv
import os, time, psutil

import ansys.mapdl.core as pymapdl

def find_process_pid(process_name):
    pid = None
    for proc in psutil.process_iter(['pid', 'name']):
        if proc.info['name'] == process_name:
            pid = proc.info['pid']
            break
    return pid

def compression(mesh,force,count):

    original_path = os.getcwd()
    relative_path = './data/'

    # managing errors from launching PyMAPDL
    try:
        start_time = time.time()
        mapdl = pymapdl.launch_mapdl(
            #jobname="Database Test",  #causes errors if turned on
            nproc=1,
            #loglevel="ERROR"
            start_timeout=60.0,
            cleanup_on_exit=True,
            port=88000,
            additional_switches='-smp'
        )
        print(f"### {count} MAPDL LAUNCH SUCCESS")
        print(f"### {count} SIMULATION STARTING")

        mesh_path = os.path.join(original_path, relative_path,str(count),mesh)
        # print("MESH PATH: "+mesh_path)
        mapdl.input(mesh_path)
        print(mapdl.cmlist())

        mapdl.prep7()
        mapdl.units("SI", mute=True)

        # eplot_path = os.path.join(original_path, relative_path,str(count),"eplot.jpg")
        # mapdl.eplot(vtk=True, show_edges=True, line_width=2, background="w",savefig="eplot_path")

        # Define a material
        mapdl.mp("EX", 2, 110e9)  # Elastic moduli in Pa (kg/(m*s**2))
        mapdl.mp("DENS", 2, 4510)  # Density in kg/m3
        mapdl.mp("NUXY", 2, 0.34)  # Poisson's Ratio
        mapdl.emodif('ALL', 'MAT', 2)

        # ALTERNATIVE: Tolerance Selection (not very good)
        # tol=0.0000108
        # base_nodes = mapdl.nsel("S", "LOC", "Z", 0-tol, 0+tol)
        # print(f"Base Nodes: {len(base_nodes)}")
        # force_nodes = mapdl.nsel("S", "LOC", "Z", 0.0045-tol, 0.0045+tol)
        # print(f"Force Nodes: {len(force_nodes)}")
        # mapdl.cm("FORCE","NODE")

        # Fix base.
        base_nodes = mapdl.components["NODESET2"].items
        print(f"Base Nodes: {len(base_nodes)}")
        mapdl.d("NODESET2", "ALL", 0)

        # Force application
        mapdl.cp(1, "UZ", "NODESET1") #reference num, degrees of freedom, which nodes
        force_nodes = mapdl.components["NODESET1"].items
        print(f"Force Nodes: {len(force_nodes)}")
        nodal_force = -force / len(force_nodes)
        mapdl.f("NODESET1", "FZ", nodal_force)

        # mapdl.nplot(plot_bc=True)

        # selecting all nodes again to solve the entire solution
        mapdl.allsel(mute=True)
        mapdl.run("/SOLU")
        mapdl.antype("STATIC")
        mapdl.solve()
        mapdl.finish()

        mapdl.post1()
        mapdl.set(1)
        
        mapdl.allsel()
        max_strain = np.nanmax(mapdl.post_processing.nodal_total_eqv_strain())
        print(max_strain)
        strain_plot_path = os.path.join(original_path, relative_path,str(count),"strain_plot.jpg")
        mapdl.post_processing.plot_nodal_total_eqv_strain(
            cpos="iso",
            lighting=False,
            background="white",
            edge_color="black",
            show_edges=True,
            n_colors=9,
            cmap="jet",
            # return_plotter=True,
            savefig=strain_plot_path
        )
        # plot.add_text('EQV Strain')
        # plot.show()

        mapdl.allsel()
        max_stress = np.nanmax(mapdl.post_processing.nodal_eqv_stress())
        print(max_stress)
        stress_plot_path = os.path.join(original_path, relative_path,str(count),"stress_plot.jpg")
        mapdl.post_processing.plot_nodal_eqv_stress(
            cpos="iso",
            lighting=False,
            background="white",
            edge_color="black",
            show_edges=True,
            n_colors=9,
            cmap="jet",
            # return_plotter=True,
            savefig=stress_plot_path
        )
        # plot.add_text('EQV Stress')
        # plot.show()

        # ALTERNATIVE LEGACY SUPPORT

        # result = mapdl.result
        # result.plot_principal_nodal_stress(
        #     0,
        #     "SEQV",
        #     lighting=False,
        #     background="white",
        #     show_edges=True,
        #     text_color="k",
        #     add_text=False,
        # )
        # nnum, stress = result.principal_nodal_stress(0)
        # von_mises = stress[:, -1]  # von-Mises stress is the right most column
        # max_stress = np.nanmax(von_mises)
        # print(max_stress)

        mapdl.finish()
        mapdl.exit()

        print(f"### {count} SIMULATION COMPLETED")
        print(f"--- %s seconds ---" % (time.time() - start_time)) #prints completed time
        return max_stress, max_strain

    except Exception as e:
        print(f"### {count} MAPDL LAUNCH ERROR ENCOUNTERED: {e}")
        print(f"--- %s seconds ---" % (time.time() - start_time)) #prints completed time

        # I have turned off the termination function to avoid interfering with any other software that is running
        # will quit() after encountering an error
        
        quit()

        # pid = find_process_pid("ANSYS.exe")

        # if pid:
        #     try:
        #         process = psutil.Process(pid)
        #         process.terminate()  # Gracefully terminate the process
        #         process.wait(timeout=5)  # Wait for the process to terminate
        #         print(f"### {count} Process {pid} terminated.")
        #     except psutil.NoSuchProcess:
        #         print(f"### {count} Process {pid} does not exist.")
        #     except psutil.AccessDenied:
        #         print(f"### {count} Permission denied to terminate process {pid}.")
        #     except psutil.TimeoutExpired:
        #         print(f"### {count} Process {pid} did not terminate in time.")
        #     finally:
        #         quit()
        # else:
        #     print(f"### {count} Process PID was not found.")

def compressionV2(mesh,force,count):

    original_path = os.getcwd()
    relative_path = './data/'
    
    log_path = os.path.join(original_path, relative_path,str(count),"pymapdl_compressionV2.txt")
    # managing errors from launching PyMAPDL
    try:
        start_time = time.time()
        mapdl = pymapdl.launch_mapdl(
            #jobname="Database Test",  #causes errors if turned on
            nproc=1,
            # loglevel="DEBUG",
            start_timeout=60.0,
            cleanup_on_exit=True,
            # log_apdl=log_path,
            additional_switches='-smp',
            port=80808,
            # exec_file=r"C:\Program Files\ANSYS Inc\ANSYS Student\v241\ansys\bin\winx64\ANSYS241.exe", # version 24.1
            version=23.2
        )
        print(f"### {count} MAPDL LAUNCH SUCCESS")
        print(f"### {count} SIMULATION STARTING")

        mapdl.clear()

        mesh_path = os.path.join(original_path, relative_path,str(count),mesh)
        print("MESH PATH: "+mesh_path)
        mapdl.input(mesh_path)
        print(mapdl.cmlist())

        mapdl.prep7()
        mapdl.units("SI", mute=True)

        # Define a material
        mapdl.mp("EX", 2, 110e9)  # Elastic moduli in Pa (kg/(m*s**2))
        mapdl.mp("DENS", 2, 4510)  # Density in kg/m3
        mapdl.mp("NUXY", 2, 0.34)  # Poisson's Ratio
        mapdl.emodif('ALL', 'MAT', 2)

        eplot_path = os.path.join(original_path, relative_path,str(count),"eplot.jpg")
        mapdl.eplot(vtk=True, show_edges=True, line_width=2, background="w",savefig=eplot_path)

        # Fix base.
        base_nodes = mapdl.components["NODESET2"].items
        print(f"Base Nodes: {len(base_nodes)}")
        mapdl.d("NODESET2", "ALL", 0)

        # Force application
        mapdl.cp(1, "UZ", "NODESET1") #reference num, degrees of freedom, which nodes
        force_nodes = mapdl.components["NODESET1"].items
        print(f"Force Nodes: {len(force_nodes)}")
        nodal_force = -force / len(force_nodes)
        mapdl.f("NODESET1", "FZ", nodal_force)

        # mapdl.nplot(plot_bc=True)

        # selecting all nodes again to solve the entire solution
        mapdl.allsel(mute=True)
        mapdl.run("/SOLU")
        mapdl.antype("STATIC")
        mapdl.solve()
        mapdl.finish()

        mapdl.post1()
        mapdl.set(1)

        mapdl.allsel()
        max_strain = np.nanmax(mapdl.post_processing.nodal_total_eqv_strain())
        print(max_strain)
        strain_plot_path = os.path.join(original_path, relative_path,str(count),"strain_plot.jpg")
        mapdl.post_processing.plot_nodal_total_eqv_strain(
            cpos="iso",
            lighting=False,
            background="white",
            edge_color="black",
            show_edges=True,
            n_colors=9,
            cmap="jet",
            font_size="20",
            # return_plotter=True,
            savefig=strain_plot_path
        )

        mapdl.allsel()
        max_stress = np.nanmax(mapdl.post_processing.nodal_eqv_stress())
        print(max_stress)
        stress_plot_path = os.path.join(original_path, relative_path,str(count),"stress_plot.jpg")
        mapdl.post_processing.plot_nodal_eqv_stress(
            cpos="iso",
            lighting=False,
            background="white",
            edge_color="black",
            show_edges=True,
            n_colors=9,
            cmap="jet",
            # return_plotter=True,
            savefig=stress_plot_path
        )

        mapdl.finish()
        #mapdl.exit()

        print(f"### {count} SIMULATION COMPLETED")
        print(f"--- %s seconds ---" % (time.time() - start_time)) #prints completed time
        return max_stress, max_strain
    except Exception as e:
        print(f"### {count} MAPDL LAUNCH ERROR ENCOUNTERED: {e}")
        print(f"--- %s seconds ---" % (time.time() - start_time)) #prints completed time
        
        try:
            mapdl.finish()
            mapdl.exit()
        except:
            pass

        # I have turned off the termination function to avoid interfering with any other software that is running
        # will quit() after encountering an error
        
        quit()

def compressionV2_24R1(mesh,force,count):

    original_path = os.getcwd()
    relative_path = './data/'
    
    log_path = os.path.join(original_path, relative_path,str(count),"pymapdl_compressionV2.txt")
    # managing errors from launching PyMAPDL
    try:
        start_time = time.time()
        mapdl = pymapdl.launch_mapdl(
            #jobname="Database Test",  #causes errors if turned on
            nproc=1,
            # loglevel="DEBUG",
            start_timeout=60.0,
            cleanup_on_exit=True,
            # log_apdl=log_path,
            additional_switches='-smp',
            port=80808,
            exec_file=r"C:\Program Files\ANSYS Inc\ANSYS Student\v241\ansys\bin\winx64\ANSYS241.exe", # version 24.1
            # version=232
        )
        print(f"### {count} MAPDL LAUNCH SUCCESS")
        print(f"### {count} SIMULATION STARTING")

        mapdl.clear()

        mesh_path = os.path.join(original_path, relative_path,str(count),mesh)
        print("MESH PATH: "+mesh_path)
        mapdl.input(mesh_path)
        print(mapdl.cmlist())

        mapdl.prep7()
        mapdl.units("SI", mute=True)

        # Define a material
        mapdl.mp("EX", 2, 110e9)  # Elastic moduli in Pa (kg/(m*s**2))
        mapdl.mp("DENS", 2, 4510)  # Density in kg/m3
        mapdl.mp("NUXY", 2, 0.34)  # Poisson's Ratio
        mapdl.emodif('ALL', 'MAT', 2)

        eplot_path = os.path.join(original_path, relative_path,str(count),"eplot.jpg")
        mapdl.eplot(vtk=True, show_edges=True, line_width=2, background="w",savefig=eplot_path)

        # Fix base.
        base_nodes = mapdl.components["NODESET2"].items
        print(f"Base Nodes: {len(base_nodes)}")
        mapdl.d("NODESET2", "ALL", 0)

        # Force application
        mapdl.cp(1, "UZ", "NODESET1") #reference num, degrees of freedom, which nodes
        force_nodes = mapdl.components["NODESET1"].items
        print(f"Force Nodes: {len(force_nodes)}")
        nodal_force = -force / len(force_nodes)
        mapdl.f("NODESET1", "FZ", nodal_force)

        # mapdl.nplot(plot_bc=True)

        # selecting all nodes again to solve the entire solution
        mapdl.allsel(mute=True)
        mapdl.run("/SOLU")
        mapdl.antype("STATIC")
        mapdl.solve()
        mapdl.finish()

        mapdl.post1()
        mapdl.set(1)

        mapdl.allsel()
        max_strain = np.nanmax(mapdl.post_processing.nodal_total_eqv_strain())
        print(max_strain)
        strain_plot_path = os.path.join(original_path, relative_path,str(count),"strain_plot.jpg")
        mapdl.post_processing.plot_nodal_total_eqv_strain(
            cpos="iso",
            lighting=False,
            background="white",
            edge_color="black",
            show_edges=True,
            n_colors=9,
            cmap="jet",
            # return_plotter=True,
            savefig=strain_plot_path
        )

        mapdl.allsel()
        max_stress = np.nanmax(mapdl.post_processing.nodal_eqv_stress())
        print(max_stress)
        stress_plot_path = os.path.join(original_path, relative_path,str(count),"stress_plot.jpg")
        mapdl.post_processing.plot_nodal_eqv_stress(
            cpos="iso",
            lighting=False,
            background="white",
            edge_color="black",
            show_edges=True,
            n_colors=9,
            cmap="jet",
            # return_plotter=True,
            savefig=stress_plot_path
        )

        mapdl.finish()
        #mapdl.exit()

        print(f"### {count} SIMULATION COMPLETED")
        print(f"--- %s seconds ---" % (time.time() - start_time)) #prints completed time
        return max_stress, max_strain
    except Exception as e:
        print(f"### {count} MAPDL LAUNCH ERROR ENCOUNTERED: {e}")
        print(f"--- %s seconds ---" % (time.time() - start_time)) #prints completed time
        
        try:
            mapdl.finish()
            mapdl.exit()
        except:
            pass

        # I have turned off the termination function to avoid interfering with any other software that is running
        # will quit() after encountering an error
        
        quit()
    
if __name__ == "__main__":
    compressionV2("tpms_mesh.cdb", 100, 0)
    # compressionV2("tpms_mesh.cdb", 100, 1)
    # compressionV2("tpms_mesh.cdb", 100, 2)
