from INPUT import*

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
# part_1 = GGG + unpaired
def part_1(recognition, start, stop):
    part_1 = "GGG" + recognition[start:stop]
    return part_1

# part_2 = paired + GUA + GUGUGU
def part_2(recognition, stop, end):
    part_2 = recognition[stop:end] + "GUA" + "GUGUGU"
    return part_2

# part_3 = RBS
def part_3(RBS):
    return RBS

# part_4 = reverse complement of paired
def part_4(trigger, stop, end):
    trigger_9 = trigger[stop:end]
    part_4 = "ACACAC" + "AUG" + trigger_9.reverse_complement()
    return(part_4)

# part_5 = linker
def part_5(linker):
    return linker

# part_6 = GFP
def part_6(GFP):
    return GFP

# full toehold construction
def toehold(trigger, start, stop, end, RBS, linker, GFP):
    toehold = part_1(trigger, start, stop) + part_2(trigger, stop, end) + part_3(RBS) + part_4(trigger, stop, end) + part_5(linker) + part_6(GFP)
    toehold = str(toehold).replace("T","U")
    return Seq(toehold)

# checking for stop condons
def check_stop(part):
    result = "GOOD There is no stop codon\n"
    i=0
    while(i<len(part)):
        if(part.find("UAA")%3==0 or part.find("UGA")%3==0 or part.find("UAG")%3==0):
            result="BAD There is a stop codon.\n"
            break
        i += 3
        part=part[i:len(part)-1]
    return result

# checking for restriction sites
def check_restriction(part):
    result = "GOOD There is no restriction site.\n"
    i=0
    while(i<len(part)):
        if(part.find("GAAUUC")!=-1 or part.find("UCUAGA")!=-1 or part.find("ACUAGU")!=-1 or part.find("CUGCAG")!=-1 or part.find("GCGGCCGC")!=-1):
            result= "BAD There is a restriction site.\n"
            break
        i+=1
        part=part[i:len(part)-1]

    return result

#========================================================================================================================
def analysis(dataframe,triggers):
    #change the paths before running
    path = 'C:\\Users\\Δημητρης\\Desktop\\Vienna\\RNA secondary structure prediction\\ViennaRNA Package'#Analysis output folder path
    viennaRNA_path = 'C:\\Users\\Δημητρης\\Desktop\\Vienna\\RNA secondary structure prediction\\ViennaRNA Package'#ViennaRNA folder path
    
    os.chdir( viennaRNA_path )#<----to run the viennaRNA .exes
    
    data = dataframe.copy()
    seq1 = []
    triggersRNA =[]
    for i in range(len(data.Seq1)):
        seq1.append(data.Seq1[i])
        
    indexes = []
    for i in range(1,len(triggers)+1):
        indexes.append(seq1.index('_'+ str(i))) 
    
     
    #print(indexes)    
    if not os.path.exists( path + '\\linker_input'):
        os.makedirs( path + '\\linker_input')

    if not os.path.exists( path + '\\linker_output'):
        os.makedirs( path + '\\linker_output')
    if  os.path.isfile( path+'\\linker_output\\rbslinker.zip'):
        os.remove(path+'\\linker_output\\rbslinker.zip')
    
            
    seq1_rbslinker = []
    
    for i in range(len(seq1)):
        ind = seq1[i].find(RBS)
        seq1_rbslinker.append(seq1[i][ind:])
    

    seq2 = []
    f = open('input2.txt','w+')
    j = 0
    for i in range(indexes[0]):
        f.write(seq1_rbslinker[i]+'\n')
        seq2.append(seq1[i])
        triggersRNA.append(trigger[0])
    while len(indexes)-j >= 2:
        if len(indexes) == 1:
            break
        else:
            for i in range(indexes[j]+1,indexes[j+1]):
                f.write(seq1_rbslinker[i]+'\n')
                seq2.append(seq1[i])
                triggersRNA.append(trigger[j])
            
        j+=1
        
       
    f.close()  
    

    
    
    command = 'RNAfold.exe>output2.txt'
    #Write the command , and use a text file as input
    process = subprocess.Popen([command ,'-p','input2.txt' ],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

    stdout , stderr = process.communicate()

    stdout,stderr

    if os.path.isfile( path + '\\linker_input\\input2.txt'):
            os.remove( path + '\\linker_input\\input2.txt')
        
    os.rename(viennaRNA_path +'\\input2.txt', path +'\\linker_input\\input2.txt')

    if os.path.isfile( path + '\\linker_output\\output2.txt'):
            os.remove( path + '\\linker_output\\output2.txt')
        
    os.rename(viennaRNA_path +'\\output2.txt', path +'\\linker_output\\output2.txt')
    

   

    with open(path+'\\linker_output\\output2.txt') as f:
        data2=f.readlines()
    
        DeltaG_RBS_Linker = []
        frequencies = []
    
        for j in range(1,len(data2),5):
            DeltaG_RBS_Linker.append(float(data2[j][-8:-2]))
        for j in range(4,len(data2),5):
            ind1 = data2[j].find(';')
            ind2 = data2[j].find('ensemble')
            frequencies.append(float(data2[j][ind2+9:ind1]))

        dat = { 'Sequences':seq2,'Frequency of the structure':frequencies,'DeltaG_RBSlinker':DeltaG_RBS_Linker}

        rbslinker_df = pd.DataFrame(dat)

        rbslinker_df_new = rbslinker_df.sort_values(by=['DeltaG_RBSlinker','Frequency of the structure'],ascending=False)

        compression_opts = dict(method='zip',
                                archive_name='rbslinker.csv')
        rbslinker_df_new.to_csv('rbslinker.zip', index=False,
                  compression=compression_opts)
        os.rename(path+'\\rbslinker.zip',path+'\\linker_output\\rbslinker.zip')

    

    'RNAfold_Toeholds'

    if not os.path.exists(path + '\\toehold_fold_input'):
        os.makedirs(path+'\\toehold_fold_input')
    if not os.path.exists(path +'\\toehold_fold_output'):
        os.makedirs(path+'\\toehold_fold_output')
    if  os.path.isfile( path+'\\toehold_fold_output\\mfe.zip'):
        os.remove(path+'\\toehold_fold_output\\mfe.zip')

        
    f = open('input1.txt','w+')
    j=0
    for i in range(indexes[0]):
        f.write(seq1[i]+'\n')
    while len(indexes)-j >= 2:
        if len(indexes) == 1:
            break
        else:
            for i in range(indexes[j]+1,indexes[j+1]):
                f.write(seq1[i]+'\n')
            
        j+=1
        
       
    f.close()  

    command = 'RNAfold.exe>output1.txt'
    #Write the command , and use a text file as input
    process = subprocess.Popen([command ,'-p','input1.txt' ],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

    stdout , stderr = process.communicate()

    stdout,stderr

    if os.path.isfile( path + '\\toehold_fold_input\\input1.txt'):
            os.remove( path + '\\toehold_fold_input\\input1.txt')
        
    os.rename(viennaRNA_path +'\\input1.txt', path +'\\toehold_fold_input\\input1.txt')

    if os.path.isfile( path + '\\toehold_fold_output\\output1.txt'):
            os.remove( path + '\\toehold_fold_output\\output1.txt')
        
    os.rename(viennaRNA_path +'\\output1.txt', path +'\\toehold_fold_output\\output1.txt')

    with open(path+'\\toehold_fold_output\\output1.txt') as f:
        data2=f.readlines()
        Sequences = []
        Energies = []
        frequencies_mfe = []
        for j in range(0,len(data2),5):
            Sequences.append(data2[j][:-1])
        for j in range(1,len(data2),5):
            Energies.append(float(data2[j][-8:-2]))
        for j in range(4,len(data2),5):
            ind1 = data2[j].find(';')
            ind2 = data2[j].find('ensemble')
            frequencies_mfe.append(float(data2[j][ind2+9:ind1]))

        dat = { 'Sequences':Sequences,'Frequency of the structure':frequencies_mfe,'MFE':Energies}

        mfe_df = pd.DataFrame(dat)

        mfe_df_new = mfe_df.sort_values(by=['MFE'],ascending=True)

        compression_opts = dict(method='zip',
                                archive_name='mfe.csv')
        mfe_df_new.to_csv('mfe.zip', index=False,
                  compression=compression_opts)
        os.rename(path+'\\mfe.zip',path+'\\toehold_fold_output\\mfe.zip')

    
    

    '''duplex'''

    #folder RNAduplex-perfect_matches
    if not os.path.exists( path + '\\duplex_input'):
        os.makedirs(path + '\\duplex_input')

    if not os.path.exists( path + '\\duplex_output'): 
        os.makedirs(path + '\\duplex_output')

    if  os.path.isfile( path+'\\duplex_output\\out.zip'):
        os.remove(path+'\\duplex_output\\out.zip')

    f = open('input3.txt','w+')

    j = 0
    for i in range(indexes[0]):
        f.write(seq1[i]+'\n'+ trigger[0] +'\n\n')
    while len(indexes)-j >= 2:
        if len(indexes) == 1:
            break
        else:
            for i in range(indexes[j]+1,indexes[j+1]):
                f.write(seq1[i]+'\n'+ trigger[j+1] +'\n\n')
            
        j+=1
 
    f.close()

    command = 'RNAduplex.exe<input3.txt>output3.txt'
    #Write the command , and use a text file as input
    process = subprocess.Popen([command , 'input3.txt' ],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

    stdout , stderr = process.communicate()

    if os.path.isfile( path + '\\duplex_input\\input3.txt'):
            os.remove( path + '\\duplex_input\\input3.txt')

    os.rename(viennaRNA_path +'\\input3.txt', path +'\\duplex_input\\input3.txt')

    if os.path.isfile( path + '\\duplex_output\\output3.txt'):
           os.remove( path + '\\duplex_output\\output3.txt')

    os.rename(viennaRNA_path +'\\output3.txt', path +'\\duplex_output\\output3.txt')
    
    with open( path +'\\duplex_output\\output3.txt','r') as f:
        data1 = f.readlines()
        
        perfect_matches = []
        ind_1 = data1[0].find('&')
        ind_2 = data1[0].find(' ')

    for i in range(len(data1)):
        ind_2=data1[i].find('&')
        ind_1=data1[i].find('(')

        counts = data1[i][ind_1:ind_2].count('(')
        perfect_matches.append((counts/len(data1[i][ind_1:ind_2-1]),counts))
          
    f =open( path + '\\duplex_output\\perfect_matches.txt','w+')
    j=0

   
    for i in range(indexes[0]):
        f.write(seq1[i] + '\t' + str(perfect_matches[i]) + '\n')#IndexError: list index out of range
        
    while len(indexes)-j >= 2:
        if len(indexes) == 1:
            break
        else:
            for i in range(indexes[j]+1,indexes[j+1]):
                f.write(seq1[i] + '\t' + str(perfect_matches[i-j-1]) + '\n')
                
                
        j+=1    
    
    
    f.close()
    with open(path+'\\duplex_output\\perfect_matches.txt','r') as f:
        a = f.readlines()
        
        Sequences = []
        nucleotides = [float(x[-4:-2]) for x in a ]
        indexes_1 = [int(x.find('(')) for x in a ]
        indexes_2 = [int(x.find(',')) for x in a]
        scores = []
        for i in range(len(a)):
            scores.append(float(a[i][indexes_1[i]+1:indexes_2[i]]))
        for i in range(len(a)):
            Sequences.append(a[i][:indexes_1[i]-1])

        dat = { 'Sequences':seq2,'Paired_nucleotides':nucleotides,'Scores':scores}

        pm_df = pd.DataFrame(dat)

        pm_df_new = pm_df.sort_values(by=['Scores','Paired_nucleotides'],ascending=False)

        compression_opts = dict(method='zip',
                            archive_name='out.csv')
        pm_df_new.to_csv('out.zip', index=False,
                  compression=compression_opts)
        os.rename(path+'\\out.zip',path+'\\duplex_output\\out.zip')
    
    '''cofold plot'''


    #folder -> RNAcofold_trigger-switch
    if not os.path.exists( path + '\\cofold_input'):
        os.makedirs(path + '\\cofold_input')


    if not os.path.exists( path + '\\cofold_output'):
        os.makedirs(path + '\\cofold_output')
    
    if not os.path.exists( path + '\\cofold_output\\cofold_plots'):
        os.makedirs( path + '\\cofold_output\\cofold_plots')

    if not os.path.exists( path +'\\Concetration_depedency_plots'):
        os.makedirs (path + '\\Concetration_depedency_plots')

    #added to keep seqX[i] the same as its used else where and should not be modified.(X=1,...,4)
    seq1_cofold_plot = []
    Binding_energy = []

    
    for i in range(len(data.Seq1)):
        seq1_cofold_plot.append(seq1[i]+GFP)
        
    
    
    f = open( 'input4.txt','w+')
    g  = open('input.txt','w+')

    
    j=0
    for i in range(indexes[0]):
        f.write(seq1_cofold_plot[i]+'&' + trigger[0] +'\n')#make fo loop for tigger list
        g.write(seq1[i]+'&' + trigger[0] +'\n')#IndexError: list index out of range
    while len(indexes)-j >= 2:
        if len(indexes) == 1:
            break
        else:
            for i in range(indexes[j]+1,indexes[j+1]):
                f.write(seq1_cofold_plot[i]+'&' + trigger[j+1] +'\n')
                g.write(seq1[i]+'&' + trigger[j+1] +'\n')
        j+=1

    f.close()
    g.close() 
    
    seqtemp1 = seq1#it seems that seqtempx and seqx are linked ,so when .remove is used it also changes seqx <-----------
    
    for i in range(1,len(triggers)+1):
        seq1_cofold_plot.remove('_%d'%i+GFP)
        seqtemp1.remove('_%d'%i)
    
    
    
    
    # INPUT DATA FOR CONCETRATIONS 
    with open('stdin','w+') as g:
        i=10**(-7)
        while(i<0.2):
            g.write(str(i)+'\t'+str(i)+'\n')
            i=1.71*i
 
    with open('input.txt','r') as g:
        a = g.readlines()
    
   
    for i in range(len(seqtemp1)):#this is the part that takes a bit long to finish
        with open('temp.seq','w+') as f,open('input.txt','r')as g:
            
            a = g.readlines()
            
            f.write(a[i])
    

        process = subprocess.Popen('RNAcofold -f stdin < temp.seq > cofold.out ',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

        stdout , stderr = process.communicate()

        stdout,stderr

    

        with open('cofold.out','r') as f:
            a=f.readlines()


        A_in = []
        B_in = []
        A = []
        B = []
        AB =[]
        AA = []
        BB = []

        for k in range(9,len(a)):
            temp = []
            for j in range(len(a[k])):
        
                if a[k][j] == '\t':
                    temp.append(j)

    
            A_in.append(np.float(a[k][:temp[0]]))
            B_in.append(np.float(a[k][temp[0]+1:temp[1]]))
            AB.append(np.float(a[k][temp[1]+1:temp[2]]))
            AA.append(np.float(a[k][temp[2]+1:temp[3]]))
            BB.append(np.float(a[k][temp[3]+1:temp[4]]))
            A.append(np.float(a[k][temp[4]+1:temp[5]]))
            B.append(np.float(a[k][temp[5]+1:-2]))


        plt.plot(A_in,A,c='red',label='Toehold')
        plt.plot(A_in,B,c='black',label='Trigger')
        plt.plot(A_in,AB,c='green',label='Toehold+Trigger')
        plt.xscale("log")
        plt.grid()
        plt.title('Concetration Depedency Plot')
        plt.ylabel('rel. concentration')
        plt.xlabel('input (mol/L) ' )
        plt.legend()
        plt.savefig('toehold_%d'%i+'.png')
        plt.close()
    
        if os.path.isfile( path + '\\Concetration_depedency_plots'+'\\toehold_%d'%i+".png"):
            os.remove( path + '\\Concetration_depedency_plots'+'\\toehold_%d'%i+".png")
        
        os.rename(viennaRNA_path +'\\toehold_%d'%i+'.png',path+'\\Concetration_depedency_plots'+'\\toehold_%d'%i+".png")
    
        with open( 'input4.txt','r') as f , open( 'my_data_%d'%i+'.txt','w') as g:
            a = f.readlines()
            g.write(a[i])

        if os.path.isfile( path + '\\rna.ps'):
            os.remove( path + '\\rna.ps')
        
        if os.path.isfile( path + '\\switch_%d'%i+'.ps'):
            os.remove( path + '\\switch_%d'%i+'.ps')

        
        command = 'RNAcofold.exe>my_outdata_%d'%i+'.txt'

        process = subprocess.Popen([command  , '-p','my_data_%d'%i+'.txt' ],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

        stdout , stderr = process.communicate() 

        stdout,stderr
        
        os.rename( viennaRNA_path + '\\rna.ps',
                   viennaRNA_path + '\\switch_%d'%i+'.ps')

        with open('my_outdata_%d'%i+'.txt') as h:
            bind_data = h.readlines()
            Binding_energy.append(float(bind_data[3][-7:]))
    
    
        
        if os.path.isfile( path + '\\cofold_output\\cofold_plots\\switch_%d'%i+'.ps'):
            os.remove( path + '\\cofold_output\\cofold_plots\\switch_%d'%i+'.ps')    
        
        os.rename( viennaRNA_path +'\\switch_%d'%i+'.ps',path+'\\cofold_output\\cofold_plots\\switch_%d'%i+'.ps')
        
        if os.path.isfile( path + '\\cofold_output\\'+'my_outdata_%d'%i+'.txt'):
            os.remove( path + '\\cofold_output\\'+'my_outdata_%d'%i+'.txt')
        
        os.rename( viennaRNA_path +'\\my_outdata_%d'%i+'.txt',path+'\\cofold_output\\'+'my_outdata_%d'%i+'.txt')
    
        if os.path.isfile( path + '\\cofold_input\\my_data_%d'%i+'.txt'):
           os.remove( path + '\\cofold_input\\my_data_%d'%i+'.txt')
        
        os.rename( viennaRNA_path +'\\my_data_%d'%i+'.txt',path+'\\cofold_input\\my_data_%d'%i+'.txt')

    if os.path.isfile( path + '\\cofold_input\\input4.txt'):
             os.remove( path + '\\cofold_input\\input4.txt')

    os.rename( viennaRNA_path +'\\input4.txt',path+'\\cofold_input\\input4.txt')

    if  os.path.isfile( path+'\\cofold_output\\cofold.zip'):
        os.remove(path+'\\cofold_output\\cofold.zip')

    
    
    dat = { 'Sequences':seq2,'DeltaG_binding':Binding_energy}

    complex_df = pd.DataFrame(dat)

    complex_df_new = complex_df.sort_values(by=['DeltaG_binding'],ascending=True)

    compression_opts = dict(method='zip',
                            archive_name='cofold.csv')
    complex_df_new.to_csv('cofold.zip', index=False,
              compression=compression_opts)
    os.rename(path+'\\cofold.zip',path+'\\cofold_output\\cofold.zip')
    
    '''plot'''
 
    #folder ->RNAplot-RNAfold_toehold
    if not os.path.exists( path + '\\fold_input'):
        os.makedirs( path + '\\fold_input')
    
    if not os.path.exists( path + '\\fold_output'):
        os.makedirs( path + '\\fold_output')  
        
    #Sequences = seq1

    for i in range(len(seq1)):
    
        f =open('input_%d'%(i+5)+'.txt','w+')
        f.write(seq1[i])
        f.close()
    
        command = 'RNAfold.exe>output_%d'%(i+5)+'.txt'
        process = subprocess.Popen([command, 'input_%d'%(i+5)+'.txt'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

        stdout , stderr = process.communicate()

        stdout,stderr
     
        if os.path.isfile( path + '\\fold_input\\input_%d'%(i+5)+'.txt'):
            os.remove( path + '\\fold_input\\input_%d'%(i+5)+'.txt')
       
        os.rename(viennaRNA_path + '\\input_%d'%(i+5)+'.txt',path+'\\fold_input\\input_%d'%(i+5)+'.txt')#<-----it gets stuck here
    
        if os.path.isfile( path + '\\fold_output\\output_%d'%(i+5)+'.txt'):
            os.remove( path + '\\fold_output\\output_%d'%(i+5)+'.txt')
        
        os.rename(viennaRNA_path + '\\output_%d'%(i+5)+'.txt', path + '\\fold_output\\output_%d'%(i+5)+'.txt')
    

        with open( path + '\\fold_output\\output_%d'%(i+5)+'.txt','r') as f ,open('output_new_%d'%(i+5)+'.txt','w') as r :
            a = f.readlines()
            index_1=a[1].index(' ')
            line_1 = a[0]
            line_2 = a[1][:index_1]
            r.write(line_1)
            r.write(line_2)
        if os.path.isfile( path + '\\rna.svg'):
            os.remove( path + '\\rna.svg')
        
        if os.path.isfile( path + '\\toehold_%d'%i+'.svg'):
            os.remove( path + '\\toehold_%d'%i+'.svg')
        
        process = subprocess.Popen('cat output_new_%d'%(i+5)+'.txt | RNAplot.exe -o svg ',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

        stdout , stderr = process.communicate()

        stdout,stderr
            
        os.rename( viennaRNA_path + '\\rna.svg',
                  viennaRNA_path + '\\toehold_%d'%i+'.svg')
    
        os.rename(viennaRNA_path+'\\toehold_%d'%i+'.svg',path+'\\fold_output\\toehold_%d'%i+'.svg')
        os.rename(viennaRNA_path+'\\output_new_%d'%(i+5)+'.txt',path+'\\fold_output\\output_new_%d'%(i+5)+'.txt')


    if not os.path.exists( path + '\\final_results'):
        os.makedirs(path + '\\final_results')
    if  os.path.isfile( path+'\\final_results\\final.zip'):
        os.remove(path+'\\final_results\\final.zip')



    
    
          

    dat = { 'Sequences':seq2,'Targets':triggersRNA,'DeltaG_binding':Binding_energy,'MFE':Energies,'Frequency_mfe':frequencies_mfe,'Paired_nucleotides':nucleotides,'Scores':scores,'Delta_G_Rbslinker':DeltaG_RBS_Linker}

    complex_df = pd.DataFrame(dat)

    complex_df_new = complex_df.sort_values(by=['DeltaG_binding','MFE'],ascending=True)

    compression_opts = dict(method='zip',
                        archive_name='final.csv')
    complex_df_new.to_csv('final.zip', index=False,
            compression=compression_opts)
    os.rename(path+'\\final.zip',path+'\\final_results\\final.zip')

    os.rename(path+'\\cofold_input',path+'\\final_results\\cofold_input')
    os.rename(path+'\\cofold_output',path+'\\final_results\\cofold_output')
    os.rename(path+'\\duplex_input',path+'\\final_results\\duplex_input')
    os.rename(path+'\\duplex_output',path+'\\final_results\\duplex_output')
    os.rename(path+'\\fold_input',path+'\\final_results\\fold_input')
    os.rename(path+'\\Concetration_depedency_plots',path+'\\final_results\\Concetration_depedency_plots')
    os.rename(path+'\\fold_output',path+'\\final_results\\fold_output')
    os.rename(path+'\\linker_input',path+'\\final_results\\linker_input')
    os.rename(path+'\\linker_output',path+'\\final_results\\linker_output')
    os.rename(path+'\\toehold_fold_input',path+'\\final_results\\toehold_fold_input')
    os.rename(path+'\\toehold_fold_output',path+'\\final_results\\toehold_fold_output')

        
    return
#=======================================================================================================================

# design of the toeholds ~ many recognition sequences for each trigger and a certain value of unpaired bases 
# based on the rule: paired + unpaired
def standard_design():
    f_1 = open("results.txt","w+")
    f_2 = open("toeholds.txt","w+")
    
    df = pd.DataFrame([], columns=['Seq1'])
    ind = 1
    for j in range(len(trigger)):
        unpaired = 9
        while(paired+unpaired<=len(trigger[j])):
            f_1.write("\n\nTRIGGER_" + str(j+1)+"\t for unpaired = "+str(unpaired) +"\n\n")

            for i in range(len(trigger[j])-paired-unpaired+1):
                trigger_j = Seq(trigger[j])
                full_toehold = toehold(trigger_j.reverse_complement(), i, i+unpaired, i+paired+unpaired, RBS, linker, GFP)
                f_2.write(str(full_toehold[0:full_toehold.find(GFP)])+"\n")
                f_1.write(str(i+1) + "\n"+check_stop(full_toehold[full_toehold.find("ACACAC"):(len(full_toehold)-1)])+check_restriction(full_toehold)+ str(full_toehold) + "\n" + str(i+1) +"\n\n")
                result_1 = check_stop(full_toehold[full_toehold.find("ACACAC"):(len(full_toehold)-1)])
                result_2 = check_restriction(full_toehold)
                if result_1 == "GOOD There is no stop codon\n" and result_2 == "GOOD There is no restriction site.\n" :
                    df_temp = pd.DataFrame([str(full_toehold[0:full_toehold.find(GFP)])], columns = ['Seq1'])
                    df = df.append(df_temp)
            unpaired+=1
        
        df_temp = pd.DataFrame(['_' + str(ind)], columns = ['Seq1'])
        df = df.append(df_temp)
        ind += 1
    f_1.close()
    f_2.close()
    
    df = df.reset_index()
    analysis(df,trigger)
    
    return 

standard_design()
