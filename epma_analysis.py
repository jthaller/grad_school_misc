# Jeremy Thaller
# Aug 19, 2020 
# Turn EPMA weight percentages into atomic compositions
# The # %% between sections exists because I've set this up 
# to run as interactive code is VS Code 


# Import the spreadsheet and split into 2 dataframes
# %%
FILE_NAME = 'AGA/IRAP08.xlsx'
SHEET = 'IRAP08'
import pandas as pd
df = pd.read_excel(io=FILE_NAME, sheet_name=SHEET, skiprows=1)
# df.head(18) 
#rows 1:15 night line
#rows 16:-1 day points
df_nline = df[0:15]
df_dpoints = df[15:]
df_dpoints.head()

# Modified from atomic_percent_converter.py
#%%
cat_wts = {"NA": 22.989769,
           "MG": 24.305,
           "SI": 28.0855,
           "AL": 26.981539,
           "P":  30.973762,
           "K":  39.0983,
           "CA": 40.078,
           "TI": 47.867,
           "FE": 55.845,
           "CR": 51.9961,
           "MN": 54.938044,
           "S":  32.065,
           "CL": 35.453}

'''
Takes a slice of the dataframe and returns an array with the mean weight
of each oxide, normalized so that the total of all weight percents = 100%
Inputs: dataframe with each column as oxide wt percent, and rows as trial num
Outputs: array of wt percents
'''
#%%
def mean_and_normalize(oxides_df):
    oxides_df = oxides_df.drop(oxides_df.columns[0:2], axis=1)
    #create new row = mean
    oxides_df.loc['mean'] = oxides_df.mean()
    #get rid of everything other than oxide weights
    oxides_df = oxides_df.loc['mean'][0:14].values
    #normalize so percentages add up to 100% total
    return [grain/oxides_df[-1] *100 for grain in oxides_df]

#%% # GRAIN 1 NIGHTTIME LINE
grain_1_ln = mean_and_normalize(df_nline)
print(grain_1_ln)

#%% # DAYTIME, SCATTER POINTS
grain_1_pts = mean_and_normalize(df_dpoints[0:5])
grain_2 = mean_and_normalize(df_dpoints[5:8])
grain_3 = mean_and_normalize(df_dpoints[8:13])
grain_4 = mean_and_normalize(df_dpoints[13:16])
grain_5 = mean_and_normalize(df_dpoints[16:20])
grain_6 = mean_and_normalize(df_dpoints[20:23])
grain_7 = mean_and_normalize(df_dpoints[23:26])
grain_8 = mean_and_normalize(df_dpoints[26:29])
grain_9 = df_dpoints[29:33]
#skip row 32 because it's pure SiO2 in grain 9.
#maybe skip row 33 too because it only has 69% total weight
grain_9.drop(grain_9.index[2], inplace=True)
grain_9 = mean_and_normalize(grain_9)
grain_10 = mean_and_normalize(df_dpoints[-1:])
#grain 10 only has 36% weight total. Maybe drop this grain
print(df_dpoints[16:20].info)
# print(grain_10)



# %%
#list of namesin the same order as above
cats = [cat for cat  in cat_wts]
# In the order above, eg [(Na2O), (MgO), etc.]
oxide_nums = [(2,1),(1,1), (1,2), (2,3), (2,5), (2,1), (1,1), (1,2), (1,1), (2,3), (1,1), (1,3), (1,0)]

'''
converts sample's weight percents to mols of each atom
inputs: cat_weight = dict {key='atom': val=amu weight}
        num_cat = int = number of cations in the oxide e.g. 1 for SiO2
        num_o = int = number of oxygens in oxide e.g. 2 for SiO2
        weight_percent = array of floats [weight percentages of each oxide measured via EPMA]
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
inputs: sample = array of weight percents for each oxide in the sample
returns: dict {key=atom: val=mols in sample}
'''
def calc_mols(sample):
    samp_mols = {"O": 0}
    for i in range(0, len(sample)-1):
        cat, o = convert(cat_wts[cats[i]], oxide_nums[i][0], oxide_nums[i][1], sample[i])
        samp_mols[cats[i]] = cat
        samp_mols["O"] += o
    return samp_mols

#run everything. These are all dicts
#keys = string of atom's name
#vals = mols of atom in the sample number
grain_1_ln_mols = calc_mols(grain_1_ln) 
grain_1_pts_mols = calc_mols(grain_1_pts)
grain_2_mols = calc_mols(grain_2)
grain_3_mols = calc_mols(grain_3)
grain_4_mols = calc_mols(grain_4)
grain_5_mols = calc_mols(grain_5)
grain_6_mols = calc_mols(grain_6)
grain_7_mols = calc_mols(grain_7)
grain_8_mols = calc_mols(grain_8)
grain_9_mols = calc_mols(grain_9)
grain_10_mols = calc_mols(grain_10)

'''
convert the mols to atoms, normalized by number of O atoms (6 or 8)
inputs: samp - a dict with keys = 'atom', and val = floats of mols
        int samp_num = which EMPA sample is being tested. 
        input num_O = number of Oxygen atoms to which you want the 
                      chemical formula normalized
returns: none. Prints table of atom names and atom numbers in chemical formula
        before and after normalization
'''
def atomic_nums(samp, samp_num=0, num_O=6):
    #first print un-normalized
    print(f"\n---Grain: {samp_num}---\nAtom\tMols")
    for atom, mol in samp.items():
        print('{}\t{}'.format(atom, round(mol,2)))
    #now normalize and print
    print(f"\n---Grain: {samp_num} Normalized---\nAtom\tMols")
    for atom, mol in samp.items():
        norm = (num_O / samp["O"])
        mol = mol*norm
        print('{}\t{}'.format(atom, round(mol,2)))

#call this to run. see above for params
# atomic_nums(grain_1_ln_mols, samp_num=1, num_O=8)
# atomic_nums(grain_1_pts_mols, 1, 8)
# atomic_nums(grain_2_mols, 2, 8)
atomic_nums(grain_10_mols, 10, 8)

# %%
# -----RESULTS-----
'''
Grain_1_ln = NaAlSi3O8 but with less Na = Albite. Small Ca inclusion
Grain_1_pts =  KAlSi3O8 = Orthoclase. But not as much K measured and a modicum of Na inclusion.
Grain_2 = KAlSi3O8 = Orthoclase. Basically the same as Grain_1_pts. I think he screwed up the labeling.
The weight percentages for grain 3 are like .01 different than those for grain 1_pts. Not an error they are the same
exact atomic percentages when normalized.

Grain_3 = Pure SiO2 = Quartz
Grain_4 = Pure SiO2 = Quartz
# GRAIN 5 data is confusing. not sure what it is. Come back.
Grain_5 = FeAlSi2O8 with .27 of Mg and .45 of K, .02 of Mn, and .04 of Cl

Grain_6 = same exact atomic percentages as above
Grain_7 = KAlSi3O8 = Orthoclase. But not as much K measured and a modicum of Na inclusion.
Grain_8 = KAlSi3O8 = Orthoclase. But 1/2 K measured and a modicum of Na inclusion.
Grain_9 = KAlSi3O8 = Orthoclase. But less Si, .4K, .15 Ca, .05 Ti, and tiny bit of Mn and Cl (.01)
Maybe redo grain 9 because it's weird and the weight percentages only added up to 69%
Grain_10 = SiO2 = Quartz. But with a modicum of Fe and Cr inclusions.
grain 10 data only had total weight percent of 36%. Not sure if this is reliable data
'''