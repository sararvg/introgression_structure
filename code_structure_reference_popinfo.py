import random
from subprocess import Popen
from shutil import copyfile
import glob
import os
import numpy as np
import re

filepath = '/opt/structure/' # path for Structure
mpath = '/home/structure/'  # main path
path_pop = 'populations/Fxx_Ax/'  # path of populations, change according to the population set you want to use
m = 0.05      # migration rate
Ngen = 10     # number of generations
Nind = 1000   # number of individuals for each generation
Nmarker = 30 # number of markers to be analysed
N1 = 70      # number of individuals from pop A
N3 = 30        # number of reference individuals from ancestral pop A
N2 = 100      # number of individuals from pop B
fst = 0.1    # Fst between the two initial populations
allele_number = 5  # number of alleles
num_individuals = range(1, N1 + N2 + N3 + 1)  # list of number of indviduals to be analysed

#1: generation 0: create a subsample with 1000 individuals of each population
#2: add to the matrix of genotypes the q value in an additional column
num_pop = 0

for dat in glob.glob(os.path.join(mpath + path_pop, '*.dat')):
    with open(dat, 'r') as population:      # open each file
        num_pop = num_pop + 1
        w = 0
        i = 0
        pop = []
        for line in population:
            w = w + 1
            line = line.replace('   ', '').replace('\t', ' ').rstrip()
            if w != 1:
                if line.startswith('1') or line.startswith('2'):
                    i = i + 1
                    if i == 1:
                        pop = line
                    else:
                        pop = np.vstack([pop, line])
    np.savetxt(dat, pop, fmt='%s', delimiter=' ')
    pop = np.genfromtxt(dat, dtype='str')       # read as matrix
    pop1 = pop[pop[:, 0] == '1', :]             # extract population 1
    pop2 = pop[pop[:, 0] == '2', :]             # extract population 2
    pop1 = np.delete(pop1, np.s_[0], axis=1)    # delete first column: number of population
    pop2 = np.delete(pop2, np.s_[0], axis=1)
    q2 = np.zeros((len(pop2), 1))               # vector with q value of pop2
    pop2 = np.hstack([pop2, q2])                # add q value to pop2
    sub1 = []
    sub1 = pop1[np.random.choice(pop1.shape[0], 1000, replace=False), :]  # subsample population 1
    q1 = np.ones((len(sub1), 1))                                          # vector with q value of sub1
    sub1 = np.hstack([sub1, q1])                                          # add q value to sub1
    subStructure3 = sub1[np.random.choice(sub1.shape[0], N3, replace=False), :]              # subsample from generation 0

    # 3: generation 1 to 100: create 1000 individuals, choosing a father and a mother from one population
    # according to the migration rate. Take a random number from 0 to 1: if the migration rate is 5%, if the number
    # is 0.95 or less father will come from population 1, if it is bigger than 0.95 he will come from population 2.
    # 4: alleles will be chosen randomly from father or mother
    # 5: add to the last column of matrix the value of q, to be obtained as the mean of q value of father and mother

    for g in range(Ngen):
        population = []
        for i in range(Nind):
            m1 = random.uniform(0, 1)
            if m1 <= (1 - m):
                parent1 = sub1[np.random.choice(sub1.shape[0], 1), :]   # random choose parent1
            else:
                parent1 = pop2[np.random.choice(pop2.shape[0], 1), :]
            parent1 = parent1[0]
            m2 = random.uniform(0, 1)
            if m2 <= (1 - m):
                parent2 = sub1[np.random.choice(sub1.shape[0], 1), :]    # random choose parent2
            else:
                parent2 = pop2[np.random.choice(pop2.shape[0], 1), :]
            parent2 = parent2[0]

            loci1 = []
            for c in range(0, len(parent1)-1):  # number of markers
                if random.randint(0, 1) == 0:
                    p1 = parent1[c][0:2]          # random choose allele from parent1
                else:
                    p1 = parent1[c][2:4]
                if random.randint(0, 1) == 0:
                    p2 = parent2[c][0:2]          # random choose allele from parent2
                else:
                    p2 = parent2[c][2:4]
                loci1.append(str(p1)+str(p2))            # add the two alleles to the list of markers
            qparent1 = float(parent1[len(parent1) - 1])  # convert to floating the q value of parent1
            qparent2 = float(parent2[len(parent2) - 1])  # convert to floating the q value of parent2
            q = (qparent1 + qparent2)/2                  # q value of offspring
            loci1.append(q)                              # add q value to list of markers
            if i == 0:
                population = loci1
            else:
                population = np.vstack([population, loci1])   # add the list of markers to the matrix of offspring
        sub1 = population



        # 7: take a subsample of individuals to analyse with Structure
        # 8: get the vector with the q values of the 100 sampled individuals to compare it with q values obtained in Structure

    for y in range(1):
        subStructure1 = population[np.random.choice(population.shape[0], N1, replace=False), :]  # subsample for Structure without replacement
        PopData1 = [1] * N1
        PopFlag1 = [0] * N1
        subPop1 = np.vstack([PopData1, PopFlag1])
        subPop1 = np.transpose(subPop1)
        subStructure2 = pop2[np.random.choice(pop2.shape[0], N2, replace=False), :]              # subsample for Structure without replacement
        PopData2 = [2] * N2
        PopFlag2 = [0] * N2
        subPop2 = np.vstack([PopData2, PopFlag2])
        subPop2 = np.transpose(subPop2)
        PopData3 = [1] * N3
        PopFlag3 = [1] * N3
        subPop3 = np.vstack([PopData3, PopFlag3])
        subPop3 = np.transpose(subPop3)
        infoPop = np.vstack([subPop1, subPop3, subPop2])
        subStructure = np.vstack([subStructure1, subStructure3, subStructure2])                   # unify subsample of the 2 populations
        qsubsample = subStructure[:, len(subStructure[0]) - 1]                                   # extract vector with q values of subsample
        subStructure = np.delete(subStructure, len(subStructure[0]) - 1, axis=1)                 # subsample with q values removed

        qsubsample = [float(x) for x in qsubsample]
        output = 'qvalues' + '_' + str(fst) + str(allele_number) +str(m) + str(N1) + str(N3) + str(N2) + str(Nmarker) + '_' + str(num_pop)
        out_name = mpath + 'real_qvalues/'
        out_name = out_name + output
        np.savetxt(out_name, qsubsample, fmt='%f', delimiter=' ')

        positions = []
        while len(positions) != Nmarker:
            a = random.randint(0, (len(subStructure[0]) - 1))
            if a not in positions:
                positions.append(a)
            else:
                continue

        subStructure = subStructure[:, positions]
        
        sampleStructure = []
        for c in range(0, len(subStructure)):
            mat = []
            for d in range(0, Nmarker):
                a = int(subStructure[c,d][0:2])
                mat.append(a)
                b = int(subStructure[c,d][2:4])
                mat.append(b)
            if c == 0:
                sampleStructure = mat
            else:
                sampleStructure = np.vstack([sampleStructure, mat])

        sampleStructure = np.hstack([infoPop, sampleStructure])

        output1 = 'subsampleStructure' + '_' + str(fst) + str(allele_number) + str(m) + str(N1) + str(N3) + str(N2) + str(Nmarker) + '_' + str(num_pop)
        out_name1 = mpath + "input/"
        out_name1 = out_name1 + output1
        np.savetxt(out_name1, sampleStructure, fmt='%i', delimiter=' ')

path = mpath + 'input'           # path for subsamples to be analyzed with Structure
path2 = mpath + 'real_qvalues'

var1 = mpath + 'mainparams'
var2 = mpath + 'extraparams'
w = 0

table_results = ['individuals', 'real_q_values', 'Structure_q_values', 'replicates']

for filename in glob.glob(os.path.join(path, 'subsampleStructure*')):
    w = w + 1
    table1 = []
    extension = re.findall(r'([0-9.]*[0-9_]*[0-9]+)', filename)                               # get number of variables and order
    copyfile(filename, mpath + 'input_Structure/input')             # copy to folder as input
    Process = Popen(filepath + 'structure.exe %s %s' % (str(var1), str(var2),), shell=True)  # run Structure
    Process.wait()
    results = open(mpath + 'output_Structure/output_f', 'rt')       # open results of Structure
    a = []

    for line in results:
        line.rstrip()
        if '(0)    1 :' in line:
            a.append(line)                # extract lines with values
        elif '(0)    2 :' in line:
            a.append(line)

    matrix_Structure = []                 # matrix with results of Structure
    for i in range(len(a)):
        if '(0)    1 :' in a[0]:
            b = []                            # empty list to be filled with all values of each individual
            columns = []                      # empty list to be filled with q values and CI of each individual
            text = a[i]                       # select line for each individual (in range (50))
            s = text.find(':  ')
            c = text.find('\n')
            o = text[s + 3:c - 1]                 # extract only part with q values and CI
            b.append(o)                       # add values to the empty list
            b = [f.replace('     ', ' ') for f in b]
            de = b[0]
            columns.append(float(de[0:5]))    # add column 1 of q values
            if i == 0:
                matrix_Structure = columns    # create matrix with q values and CI
            else:
                matrix_Structure = np.vstack([matrix_Structure, columns])   # add to the matrix all q values and CI
        elif '(0)    2 :' in a[0]:
            b = []                            # empty list to be filled with all values of each individual
            columns = []                      # empty list to be filled with q values and CI of each individual
            text = a[i]                       # select line for each individual (in range (50))
            s = text.find(':  ')
            c = text.find('\n')
            o = text[s + 3:c - 1]                 # extract only part with q values and CI
            b.append(o)                       # add values to the empty list
            b = [f.replace('     ', ' ') for f in b]
            de = b[0]
            columns.append(1.000 - float(de[0:5]))    # add column 1 of q values
            if i == 0:
                matrix_Structure = columns    # create matrix with q values and CI
            else:
                matrix_Structure = np.vstack([matrix_Structure, columns])

    replicates = [w] * (N1 + N3 + N2)
    individual_number_list = [N1 + N3 + N2] * (N1 + N3 + N2)

    for name in glob.glob(os.path.join(path2, 'qvalues*')):
        extension2 = re.findall(r'([0-9.]*[0-9_]*[0-9]+)', name)
        if extension2 == extension:
            real_q = np.genfromtxt(name, dtype='str')
            break
        else:
            continue

    table1 = np.column_stack((num_individuals, real_q, matrix_Structure[:, 0], replicates))

    table_results = np.vstack([table_results, table1])

output2 = 'Fst' + str(fst) + '_NAl' + str(allele_number) + '_NMar' + str(Nmarker) + '_M' + str(m) + '_N' + str(N1 + N3 + N2) + '_ref_popinfo_' + str(N3)
out_name2 = mpath + "results/"
out_name2 = out_name2 + output2
np.savetxt(out_name2, table_results, fmt='%s', delimiter=' ')
