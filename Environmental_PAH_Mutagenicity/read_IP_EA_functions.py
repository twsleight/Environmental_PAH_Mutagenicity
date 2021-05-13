

import os

path = r'C:\ResearchWorkingDirectory\MutagenicityDataAug\lastfew_verts'
path = r'C:\ResearchWorkingDirectory\echa_data_gas'
filename = 'echa_2_gas.out'




def read_neut_energy(dir_path, filename):
        
    '''reads the energy from a gas phase calculation
    
    Optional: prioritized descriptors based on a list. 

    Arguments: directory path and filename
    
    Returns: energy value in Hartrees
    '''   

    with open(os.path.join(dir_path, filename)) as searchfile:    
        #************************  

        eng_val = 'not foundc'
        thermal_correction = 'not foundt'
        for line in searchfile:
            #the normal termination flag should turn to 1
            if 'Normal termination of Gaussian 09' in line:
                terminated_normally = True
            if 'slurmstepd: error:' in line:
                raise Exception('slurmstepd: error in ',filename)
            #***************************************************************
     
            #electronic energy for gas phase
            left,sep,right=line.partition('Thermal correction to Gibbs Free Energy=')
            if sep:

                t1 = line. split()
                thermal_correction = float(t1[-1])
#               #if we're unable to find  this  there's a problem     
          
            left,sep,right=line.partition('Sum of electronic and thermal Free Energies=') 
            if sep:

                c1 = line. split()
                
                #This subraction was the only way to get this data out of the raw text file
                eng_val = float(c1[-1])-thermal_correction  #this is the final data
#                    print(c)

        #error checking
        if not terminated_normally:
            raise Exception('Normal termination not detected',filename)

    return eng_val

def read_posneg_energy(dir_path, filename):     
    '''reads the energy from a vertical calculation
    
    Optional: prioritized descriptors based on a list. 

    Arguments: directory path and filename
    
    Returns: energy value in Hartrees
    '''   
    
    #read the vertical energy value
    with open(os.path.join(dir_path, filename)) as searchfile:  
        
        eng_val = 'not found'
        scf_done = 'not found'
        
        for line in searchfile:  
            
            #electronic energy for gas phase
            left,sep,right=line.partition('SCF Done:  E(UB3LYP) =')
            if sep:
                t1 = line. split()
                eng_val = float(t1[4])
                scf_done = 'good'

        #catch any computational errors
        if scf_done == 'not found':
            
            print(filename, 'has computation issues')
    
    return eng_val




