## Jeremy Thaller
## 2020.06.19
## Convert weight percentages from EPMA results into molecular formulas
import numpy as np

#dict{key = name, val = weight in amu}
cat_wts = {"NA": 22.989769,
           "MG": 24.305,
           "SI": 28.0855,
           "AL": 26.981539,
           "K":  39.0983,
           "CA": 40.078,
           "TI": 47.867,
           "FE": 55.845,
           "MN": 54.938044}

#list of namesin the same order as above
cats = ["NA", "MG", "SI", "AL", "K", "CA", "TI", "FE", "MN"]

# In the order above, eg [(Na2O), (MgO), etc.]
oxide_nums = [(2,1),(1,1), (1,2), (2,3), (2,1), (1,1), (1,2), (1,1), (1,1)]

# weight percentages from each mystery sample from the EMPA results       
sample_1 = [0.35,	15.57,	51.37,	2.37,	0.02,	20.47,	0.67,	9.16,	0.24]
sample_2 = [0.05,	24.03,	52.86,	1.03,	0.01,	2.01,	0.37,	19.30,	0.51]
sample_3 = [0.50,	2.79,	6.06,	3.46,	0.37,	0.44,	14.32,	64.65,	0.37]
sample_4 = [5.03,	0.11,	54.68,	28.10,	0.35,	11.25,	0.07,	0.82,	0.01]
sample_5 = [4.26,	2.05,	62.46,	15.36,	2.80,	4.42,	1.29,	6.20,	0.07]

'''
converts sample's weight percents to mols of each atom
inputs: cat_weight = dict {key='atom': val=amu weight}
        num_cat = int = number of cations in the oxide e.g. 1 for SiO2
        num_o = int = number of oxygens in oxide e.g. 2 for SiO2
        weight_percent = array of floats [weight percetages of each oxide measured via EPMA]
returns: dict {key=atom: val=mols in sample}
'''
def convert(cat_weight, num_cat, num_o, weight_percent):
    mol_wt_oxide = cat_weight * num_cat
    mol_oxide = weight_percent / mol_wt_oxide
    mol_cat = mol_oxide * num_cat
    mol_oxygen = mol_oxide * num_o
    return mol_cat, mol_oxygen

'''
This basically just calls the function above and takes care of the fact that O is in all the oxides
inputs: sample = arraw of weight percents for each oxide in the sample
returns: dict {key=atom: val=mols in sample}
'''
def calc_mols(sample):
    samp_mols = {"O": 0}
    for i in range(0, len(sample)-1):
        cat, o = convert(cat_wts[cats[i]], oxide_nums[i][0], oxide_nums[i][1], sample[i])
        # cat, o = convert(cat_wts[cats[i]], num_cats[i], num_o[i], sample[i])
        samp_mols[cats[i]] = cat
        samp_mols["O"] += o
    return samp_mols

#run everything. These are all dicts
#keys = string of atom's name
#vals = mols of atom in the sample number
samp_1 = calc_mols(sample_1)
samp_2 = calc_mols(sample_2)
samp_3 = calc_mols(sample_3)
samp_4 = calc_mols(sample_4)
samp_5 = calc_mols(sample_5)

'''
convert the mols to atoms, normalized by number of O atoms (6 or 8)
inputs: samp - a dict with keys = 'atom', and val = floats of mols
        int samp_num = which EMPA sample is being tested. 
        input num_O = number of Oxygen atoms to which you want the 
                      chemical formula normalized
returns: none. Prints table of atom names and atom numbers in chemical formula
        before and after normalization
'''
def atomic_nums(samp,samp_num=0, num_O=6):
    #first print un-normalized
    print(f"\n---Sample {samp_num}---\nAtom\tMols")
    for atom, mol in samp.items():
        print('{}\t{}'.format(atom, round(mol,2)))
    #now normalize and print
    print(f"\n---Sample {samp_num} Normalized---\nAtom\tMols")
    for atom, mol in samp.items():
        norm = (1 / samp["O"])*num_O
        mol = mol*norm
        print('{}\t{}'.format(atom, round(mol,2)))

#call this to run. see above for params
atomic_nums(samp_4, 4, 8)


