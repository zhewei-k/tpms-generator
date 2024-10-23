"""
NTPMS v1.0.0
------------
NTPMS is a module that assists in TPMS generation via nTopCL. Note that porosity has been turned off.
Requires the following files to be in the same folder directory
- input_template.json
- tpms.ntop
"""

#Imports
import os, subprocess, json, shutil, time, numpy as np
from volume_calculator import STLUtils


#Assuming this script, ntop file, and json files will be in the same folder
Current_Directory = os.path.dirname(os.path.realpath('__file__')) 
exePath = r"ntopcl"  #nTopCL path, raw string
nTopFilePath = r"tpms_no_porosity_dynamic.ntop"   #nTop notebook file name, raw string
input_file_name = "input_tpms.json"      #JSON input file name to be saved as
output_file_name = "output_tpms.json"       #JSON output file name to be saved as

#Input variables in JSON structure
with open('input_template.json','r') as f:
    Inputs_JSON = json.load(f)

def calculatePorosity(cellsize, count):
    """Calculates porosity

    Args:
        cellSize (float): size of unit cell (mm)
        count (int): folder reference number eg ./data/count

    Returns:
        porosity (float): volume of TPMS geometry (mm^3) / volume of bounding box (mm^3)
        mesh_density (float): number of triangles (#) / surface area (mm^2)
    """
    vcal = STLUtils() #calculate STL volume to improve speed
    elements = vcal.loadSTL('./data/'+str(count)+'/'+'tpms_mesh.stl')
    tpms_volume = vcal.calculateVolume("cm") * 1000     #volume cm^3 -> mm^3
    tpms_area = vcal.surf_area()                        #surface area cm -> mm^3
    mesh_density = elements / tpms_area                 #elements per mm^3
    bounding_box_volume = cellsize ** 3
    porosity = (bounding_box_volume - tpms_volume) / bounding_box_volume

    return porosity, mesh_density

def elementSize(cellsize, thickness, elementmultiplier):
    """Calculate element size for meshing at ideal resolution for the range of 0.1mm - 2mm TPMS thickness.

    Args:
        thickness (float): thickness of tpms unit cell (mm)        

    Returns:
        element_size (float): meshing element size
    """

    # trying to replicate a ramp followed by two fixed heights
    # The type of element sizes we want to see were first done on a 5MM unit cell, and everything is scaled from that reference.
    # element_size
    # ^      ____
    # |  ____
    # | /
    # |------------> thickness

    lower_thickness_threshold = 0.1
    upper_thickness_threshold = 0.4
    cellsize_ref = 5
    scale = cellsize_ref / cellsize
    scaled_thickness = thickness * scale    # scale to match 5MM unit cell
    
    if scaled_thickness <= lower_thickness_threshold:
        element_size = thickness * 1.1                            
    elif lower_thickness_threshold < scaled_thickness < upper_thickness_threshold:
        element_size = np.round(0.12 / scale,3) # larger element size if cell size is large
    elif scaled_thickness >= upper_thickness_threshold:
        element_size = 0.22 / scale  # larger element size if cell size is large
    else:
        raise Exception("Element size has not been successfully generated")
    return element_size * elementmultiplier

def meshList(latticetype:list,cellsize:list,thickness:list):
    """Takes in input parameters and runs nTop to output porosity and corresponding STL mesh into independent folders specified by count. Returns a summary as a CSV file.

    Args:
        latticeType (list): 0-Gyroid, 1-Schwarz, 2-Diamond, 3-Lidinoid, 4-SplitP, 5-Neovius
        cellSize (list): list of size of unit cells (mm)
        thickness (list): list of thickness of tpms unit cells (mm)
    """    
    start_time = time.time()

    #CSV File to store porosity values
    CSV_file = "tpms_meshList.csv"
    Header = ['Count', 'Lattice Type', 'Cell Size', 'Thickness', 'Porosity', 'Mesh Density']
    output_file = open(CSV_file, mode='w')
    output_file.write(",".join(Header)+'\n')
    output_file.close()

    if os.path.isdir('./data/'):    #checks whether "data" folder is present,
        shutil.rmtree('./data/')    #if yes, deletes "data" folder
    os.mkdir('./data/')             #creates folder

    count = 0
    #looping through all configurations
    for i in latticetype:
        for j in cellsize:
            for k in thickness:
                os.mkdir('./data/'+str(count))
                Inputs_JSON['inputs'][0]['value'] = i                           #looping latticeType
                Inputs_JSON['inputs'][1]['value'] = j                           #looping cellSize
                Inputs_JSON['inputs'][2]['value'] = k                           #looping thickness
                Inputs_JSON['inputs'][3]['value'] = elementSize(j,k)              #calculate element size
                Inputs_JSON['inputs'][4]['value'] = str(count)                  #ref for file location

                #defining input and output JSON paths
                input_path = './data/'+str(count)+'/'+input_file_name           #path for input JSON for that iteration
                output_path = './data/'+str(count)+'/'+output_file_name         #path for output JSON for that iteration

                #appends nTopCL arguments as a list for concatenation later
                Arguments = [exePath]               #nTopCL path
                Arguments.append("-v2")             #2nd level error report
                Arguments.append("-j")              #json input argument
                Arguments.append(input_path)        #json path
                Arguments.append("-o")              #output argument
                Arguments.append(output_path)       #output json path
                Arguments.append(nTopFilePath)      #.ntop notebook file path

                #creates input JSON file
                with open(input_path, 'w') as outfile:
                    json.dump(Inputs_JSON, outfile, indent=4)    

                #runs nTopCL with the arguments defined and saves the error output
                print(" ".join(Arguments))
                output,error = subprocess.Popen(Arguments,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

                #print all return messages
                print(output.decode("utf-8"))

                #processes output and updates the CSV file
                # with open(output_path,'r') as f:
                #     data = json.load(f)             #load json file
                #     porosity = data[0]["value"]["val"]

                porosity, mesh_density = calculatePorosity(j, count)

                summary = [count, i, j, k, porosity, mesh_density]        #attach to summary
                output_file = open(CSV_file,mode='a')       #append to csv
                output_file.write(",".join([str(i) for i in summary])+'\n')     #write to csv by inputting a piece of data and then creating a blank line
                output_file.close()                                             #close csv
                print(summary)
                count += 1

    print("--- %s seconds ---" % (time.time() - start_time)) #prints completed time

def mesh(latticetype:int,cellsize:float,thickness:float,elementmultiplier:float,count:int):
    """Takes in input parameters and runs nTop to output porosity and corresponding STL mesh into independent folders specified by count. Returns a summary of the input values as a list, including porosity.

    Args:
        latticeType (int): 0-Gyroid, 1-Schwarz, 2-Diamond, 3-Lidinoid, 4-SplitP, 5-Neovius
        cellSize (float): size of unit cell (mm)
        thickness (float): thickness of tpms unit cell (mm)
        count (int): creates a folder eg ./data/count and stores files

    Returns:
        summary (list): returns the count, latticeType, cellSize, thickness, and porosity.
    """
    start_time = time.time()

    #for running as an individual function
    if not(os.path.isdir('./data/')):    #checks whether "data" folder is present,
        os.mkdir('./data/')              #if not present, creates folder
    folder = os.path.join('./data/',str(count))
    if os.path.isdir(folder):    #checks whether "data" folder is present,
        shutil.rmtree(folder)    #if yes, deletes "data/count" folder
    os.mkdir('./data/'+str(count))

    Inputs_JSON['inputs'][0]['value'] = latticetype
    Inputs_JSON['inputs'][1]['value'] = cellsize
    Inputs_JSON['inputs'][2]['value'] = thickness
    Inputs_JSON['inputs'][3]['value'] = element_size = elementSize(cellsize, thickness, elementmultiplier)
    Inputs_JSON['inputs'][4]['value'] = str(count)                  #ref for file location

    #defining input and output JSON paths
    input_path = './data/'+str(count)+'/'+input_file_name           #path for input JSON for that iteration
    output_path = './data/'+str(count)+'/'+output_file_name         #path for output JSON for that iteration

    #appends nTopCL arguments as a list for concatenation later
    Arguments = [exePath]               #nTopCL path
    Arguments.append("-v2")             #2nd level error report
    Arguments.append("-j")              #json input argument
    Arguments.append(input_path)        #json path
    Arguments.append("-o")              #output argument
    Arguments.append(output_path)       #output json path
    Arguments.append(nTopFilePath)      #.ntop notebook file path

    print(f"### {count} TPMS Generation Started")
    print(f"### {count} Element size is {element_size}")

    with open(input_path, 'w') as outfile:          #creates input JSON file
        json.dump(Inputs_JSON, outfile, indent=4)    

    print(" ".join(Arguments))                      #runs nTopCL with the arguments defined and saves the error output
    output,error = subprocess.Popen(Arguments,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    print(output.decode("utf-8"))                   #print all return messages

    # nTop with POROSITY
    # #from output, reads porosity to update the CSV file
    # with open(output_path,'r') as f:
    #     data = json.load(f)             #load json file
    #     porosity = data[0]["value"]["val"]

    porosity, mesh_density = calculatePorosity(cellsize,count)
    summary = [count, latticetype, cellsize, thickness, porosity, element_size]        #attach to summary

    print(f"### {count} TPMS Generation Finished")
    print("--- %s seconds ---" % (time.time() - start_time))  

    return summary

if __name__ == "__main__":
    print("MAIN IS RUNNING")
    #print(mesh(0,5,1,1))
    meshList([0],[5],[0.3,0.4])

