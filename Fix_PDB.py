# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:33:17 2017


@author: amirbitran
"""
import numpy as np
import copy as cp





def fix_PDB(filepath, fixed_filepath, aa_col=5, atom_col = 1):
    """   
    The point of this script is to fix the amino acid numbering in a PDB file
    Basically loops through the PDB file and makes sure amino acid numbers are always in consecutive order and start with 1
    Be sure to check the variable aa_col, whcich indicates which column in the PDB file has the AA number (in zero-indexing)
    Also changes atom numbering if necessary
    """    
    f = open(filepath, 'r+') #name of file we want to fix
    
    f_fixed=open(fixed_filepath, 'w')  #how we want to name file once it's fixed
    
    
    n=1  #indices amino acid
    lines=[]
    l=1 #a counter for lines that start with ATOM
    for line in f.readlines():
        #print(line) 
        lines.append(line)
        line=line.rstrip('\n')
        
        if line[0:3]!='END' and line[0:4]=='ATOM':
            entries=line.split()
            aa_number=entries[aa_col]
            atom_number = entries[atom_col]
            
            if l>1:  
                if aa_number!=prev_aa_number:#amino acid has chnaged
                    n+=1
            
            prev_aa_number=cp.deepcopy(aa_number)
            
            if aa_number!=n:  #Then we need to change the line
                length_of_existing_number=len('{}'.format(aa_number))
                length_of_number=len('{}'.format(n))
                
                line_array=[c for c in line]
    
                entry_positions=[c for c in range(len(line)) if (c==0) or (line[c-1]==' ' and line[c]!=' ' )]  #where are the different entries in this line?
                aa_number_pos = entry_positions[aa_col]  #at what position in this line does the amino acid number  start?
        
    
                for c in np.arange(aa_number_pos, aa_number_pos + length_of_existing_number, 1):#clear any number that was there before
                    line_array[c]=''        
                
                line_array[aa_number_pos:aa_number_pos+length_of_number]= '{}'.format(n) 
                #line_array[25-length_of_number+1:26] = '{}'.format(n) 
                fixed_line=''.join(line_array)+'\n'     
                
            
            if atom_number !=l: #Fix atom number if it is wrong
                print(l)
                length_of_existing_number=len('{}'.format(atom_number))
                length_of_number=len('{}'.format(l))
                
                line_array=[c for c in fixed_line]
    
                entry_positions=[c for c in range(len(fixed_line)) if (c==0) or (fixed_line[c-1]==' ' and fixed_line[c]!=' ' )]  #where are the different entries in this line?
                atom_number_pos = entry_positions[atom_col]+ length_of_existing_number-1  #at what position in this line does the atom  number  END? Note numbers should be right aligned
        
    
                for c in np.arange(atom_number_pos - length_of_existing_number+1, atom_number_pos+1 , 1):#clear any number that was there before
                    line_array[c]=' '        
                
                line_array[atom_number_pos - len('{}'.format(l))+1:atom_number_pos+1]= '{}'.format(l)  #make sure numbers are right alined
                #line_array[25-length_of_number+1:26] = '{}'.format(n) 
                fixed_line=''.join(line_array)
                
            f_fixed.write(fixed_line)
                
            l+=1
                #print(fixed_line)
    
    
    f_fixed.write('END')
    f.close()    
    f_fixed.close()
    
    
def Shift_PDB_numbering(filepath, fixed_filepath, shift, starting_AA, aa_col=5):
    """
    Opens the PDB file in filepath, and shifts all amino acid numbers by a value shift (can be positive or negative)
    for all amino acids whose original number is greater than or equal to starting_AA
    
    Saves as fixed_filepath
    """
    f = open(filepath, 'r+') #name of file we want to fix
    
    f_fixed=open(fixed_filepath, 'w')  #how we want to name file once it's fixed
    
    lines=[]    

    for line in f.readlines():
        #print(line) 
        lines.append(line)
        line=line.rstrip('\n')
        
        if line[0:3]!='END' and line[0:4]=='ATOM':
            entries=line.split()
            aa_number=entries[aa_col]
            
            
            if int(aa_number)>=starting_AA:  #Then we need to change the line
                new_number = int(aa_number) + shift
            else:
                new_number = int(aa_number)
                
                
            length_of_existing_number=len('{}'.format(aa_number))
            length_of_number=len('{}'.format(new_number))
            
            line_array=[c for c in line]

            entry_positions=[c for c in range(len(line)) if (c==0) or (line[c-1]==' ' and line[c]!=' ' )]  #where are the different entries in this line?
            aa_number_pos = entry_positions[aa_col]  #at what position in this line does the amino acid number  start?
        
            
            for c in np.arange(aa_number_pos, aa_number_pos + length_of_existing_number, 1):#clear any number that was there before
                line_array[c]='' 
    
            line_array[aa_number_pos:aa_number_pos+length_of_number] = '{}'.format(new_number) 
            fixed_line=''.join(line_array)+'\n'     
            f_fixed.write(fixed_line)

    
    f_fixed.write('END')
    f.close()    
    f_fixed.close()