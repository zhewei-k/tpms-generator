'''
------------
LOGGER
------------
LOGGER keeps track of time and outputting data into a spreadsheet. 
'''

from functools import wraps
import os, shutil, time, csv
from datetime import datetime

dataset = {
    'Count':None,        #int
    'ExecutionTime':None,#str
    'ElementNum':None,   #int
    'NodeNum':None,      #int
    'LatticeType':None,  #int
    'CellSize':None,     #float
    'Thickness':None,    #float
    'Porosity':None,     #float
    'ElementSize':None,  #float
    'InternalStress':None,       #float
    'Strain':None,       #float
    'Force':None,        #float
    'Stiffness':None,    #float
}

def main():
    folder()
    print("directory created at "+ str(folder.directory))
    spreadsheet("export_test.csv")

# Defines a decorator to measure execution time
def timer(func):
    """
    Timer to record how long it takes to execute a function. Works as a decorator.
    Prints the time required to completed the function.

    Args:
        func (function): any function, it may be appropriate to add a delay if it is executed instantly.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        dataset['ExecutionTime'] = execution_time  #records the time to execute the function
        print(f"------- {func.__name__} executed in {execution_time:.4f} seconds -------\n")
        return result
    return wrapper
    

def folder(custom_directory=""):
    """Creates a folder in the data folder.
    Args:
        None
    
    Returns:
        directory (str): 
    """
    if custom_directory:
        directory = custom_directory
        print(f"### Directory reused at: {directory}")
    else:
        current_time = datetime.now().strftime("%Y_%m_%d_%H%M")
        cwd = os.getcwd()
        directory = os.path.join(cwd,"data",current_time)
        print(f"### Directory has been created at: {directory}")

        if os.path.isdir(directory):    #checks whether "data" folder is present,
            shutil.rmtree(directory)    #if yes, deletes "data" folder
        os.mkdir(directory)             #creates folder if data folder for storing data
    
    setattr(folder,'directory',directory)
    return directory


def variable(array, var):
    array.append(var)
    print(f"Array has been updated: {array}")
    

def spreadsheet(data_dict, filename='log.csv', mode='a'):
    """Logs data from a dictionary to CSV file of choice.

    Args:
        data_dict (var): Dictionary with variable names as keys and values.
        filename (str): Name of the CSV file.
        mode (str): Mode for opening the file ('a' for append, 'w' for overwrite).
    """
    abs_directory = os.path.join(folder.directory,filename)
    with open(abs_directory, mode, newline='') as csvfile:
        fieldnames = list(data_dict.keys()) #compiles dictionary keys as a list
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames) #writes header
        
        # Write header if file is new
        if csvfile.tell() == 0:
            writer.writeheader()
        
        # Write data
        writer.writerow(data_dict)

if __name__ == "__main__":
    main()