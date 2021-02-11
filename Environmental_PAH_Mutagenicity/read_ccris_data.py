# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:27:56 2019

@author: twsle
"""
#acknowledgements
#https://stackoverflow.com/a/54932071/10226776 for the CIR convert function
import pandas as pd
import ssl


# This restores the same behavior as before.
context = ssl._create_unverified_context()
# urllib.urlopen("https://no-valid-cert", context=context)

from urllib.request import urlopen

def CIRconvert(ids):
    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + ids + '/smiles'
        ans = urlopen(url, context = context).read().decode('utf8')
        return ans
    except:
        return 'Did not work'

#code for reading in the xml file

import xml.etree.ElementTree as etree
from Environmental_PAH_Mutagenicity.read_ccris_data import CIRconvert

def convert_xml_xlsx(filename):
    '''function for reading the CCRIS xml file    

    Args: filename (if i working directory) or filepath to xml file
    
    Returns: df_rows, dataframe containing the mutagenicity data
    '''   
    
    tree = etree.parse(filename)
    root = tree.getroot()

    df_rows = pd.DataFrame()
    for child in root:
        name = child.find("NameOfSubstance").text
        cas_num = child.find( "CASRegistryNumber").text
        smiles = CIRconvert(cas_num)
        ips = child.findall("mstu")

        #if ips is not found there is no data in this child
        if ips == []:

            continue
        else:
            for i in ips:

                if (i.find("matvm") is not None):
                    acts9 = i.find("matvm").text
                    if acts9 == 'NONE':

                        if (i.find("rsltm") is not None):
                            resultMut = i.find("rsltm").text

                            if (i.find('indcm') is not None):
                                resultStrain = i.find("indcm").text

                                if (i.find("tsstm") is not None):
                                    testMethod = i.find("tsstm").text
                                    #allstuff.append(acts9 + resultMut + testMethod)

                                    #build the dataframe right here
                                    tempRow = pd.DataFrame([{'CAS': cas_num,
                                        'SMILES':smiles, 'name':name,
                                        'method':acts9, 'Strain':resultStrain,
                                        'result': resultMut}])
                                    df_rows = df_rows.append(tempRow)
                    else:
                        continue

    return df_rows
