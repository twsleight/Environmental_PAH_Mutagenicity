#gethlenergy will take a filepath and a file name
#and return the HOMO and LUMO energies.

import os

def gethlenergy(dir_path, filename):

    #initilalize the error codes
    terminated_normally = False

    #initialize the variables to recieve the data
    a = 0
    b = 0
    #initialize a holding variable
    prevLine = []

    with open(os.path.join(dir_path, filename)) as searchfile:
        for line in searchfile:
            #error checking
            #the normal termination flag should turn to 1
            if 'Normal termination of Gaussian 09' in line:
                terminated_normally = True
            if 'slurmstepd: error:' in line:
                raise Exception('slurmstepd: error in ',filename)

            #The Alpha  occ. and Alpha virt lines should
            #be found right next to each other
            if 'Alpha virt. eigenvalues --' in line and 'Alpha  occ. eigenvalues --' in prevLine:
                left,sep,right=line.partition('Alpha virt. eigenvalues --')

                #This is the LUMO values
                b = float((right[:56]).split()[0])

                left,sep,right=prevLine.partition('Alpha  occ. eigenvalues --')

                #This is the HOMO value
                a = float((right[:56]).split()[-1])

            prevLine = line
        #end of the searchfile loop
        if not terminated_normally:
            raise Exception('Normal termination not detected',filename)

    return [a,b]
