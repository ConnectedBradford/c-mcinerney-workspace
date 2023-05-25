# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Helper functions
#
# This notebook collates all the functions that help the other notebooks do their thing, without clogging the other notebooks with function definitions.

# ## Imports

# + jupyter={"source_hidden": true} tags=[]
import nest_asyncio
#nest_asyncio.apply()
import numpy
import pandas
pandas.set_option('display.max_colwidth', None)
from datetime import date, datetime
import itertools
import scipy.stats
import sklearn.metrics
import math
import os
import matplotlib.pyplot
from google.cloud import bigquery, exceptions
from IPython.display import display, Markdown, Latex
from IPython import get_ipython
from tqdm import tqdm
import pyarrow.parquet
import pathlib
import timeit
import re
import rpy2.ipython
from multiprocessing.pool import Pool


# -

# ## Functions

# ### `entropy_output()`: Compute and present the entropy

# + jupyter={"source_hidden": true} tags=[]
# A function to compute and present the entropy of a column in a pandas.Dataframe.
#
# ARGUMENTS
# 1. my_df_column:            A column from a pandas.Dataframe containing the variable
#                             for which entropy needs calculating.
#
# RETURNS
# 1. entropy_caseness_scaled: Entropy scaled to the theoretical maximum for a binary variable.
# 2. entropy_caseness:        Entropy in nats.
#

def entropy_output(my_df_column):
    entropy_caseness = scipy.stats.entropy(my_df_column.value_counts(), base = math.e)
    entropy_caseness_scaled = round(entropy_caseness / math.log(2, math.e) * 100, 3)
    if entropy_caseness < 0.001:
        print("\t Caseness variable entropy < 0.001 nats")
    else:
        print("\t Caseness variable entropy =", round(entropy_caseness, 3), "nats")
    if entropy_caseness < 0.001:
        print("\t Caseness variable scaled entropy < 0.001 %")
    else:
        print("\t Caseness variable scaled entropy =", entropy_caseness_scaled, "%")
    
    return entropy_caseness, entropy_caseness_scaled


# -

# ### `hitrate_output()`: Compute and present the hit rate - a.k.a. (mis)classification error - of a caseness variable

# + jupyter={"source_hidden": true} tags=[]
# A function to compute and present the hitrate - a.k.a. (mis)classification
# error - of a caseness variable.
#
# ARGUMENTS
# 1. my_caseness_variable:  A column from a pandas.Dataframe containing patients'
#                           caseness values.
#
# RETURNS
# 1. hitRate_none:          The hitrate assuming no one demonstrates the caseness.
# 2. hitRate_all:           The hitrate assuming everyone demonstrates the caseness.
#

def hitrate_output(my_caseness_variable):
    numerator = (my_caseness_variable != 0).sum()
    denominator = len(my_caseness_variable)
    hitRate_all = (numerator / denominator) * 100
    hitRate_none = 100 - hitRate_all
    Odds_noYes = hitRate_none / hitRate_all
    if hitRate_all < 0.001:
        print("\t Hit rate (all) < 0.001 %")
        print("\t Hit rate (none) \u2248 100 %")
        print("\t Odds (No CMHD : CMHD) \u2248 infitely-times less likely to have CMHD than to have it.")
    elif hitRate_all:
        print("\t Hit rate (all) =", round(hitRate_all, 3), "%")
        print("\t Hit rate (none) =", round(hitRate_none, 3), "%")
        print("\t Odds (No CMHD : CMHD) =", f'{int(Odds_noYes):,}', "-times less likely to have CMHD than to have it.")
    else:
        print("\t Hit rate (all) =", round(hitRate_all, 3), "%")
        print("\t Hit rate (none) =", round(hitRate_none, 3), "%")
        print("\t Odds (No CMHD : CMHD) =", f'{int(Odds_noYes):,}', "-times less likely to have CMHD than to have it.")
    
    return hitRate_none, hitRate_all


# -

# ### `boundaryfilter()`: Filter out feature sets that are not within the prevalence bounds

# + jupyter={"source_hidden": true} tags=[]
# A function to filter out feature sets that are not within the prevlance bounds.
#
# The function counts the non-zero elements of each feature set can compares it 
# to the minimum and maximum count criteria.
#
# ARGUMENTS
# 1. my_featureSet_array:       An n-by-fs pandas.Datqframe of n patients represented
#                               by rows and fs features represented by columns.
# 2. verboase:                  An optional argument to indicate whether the user 
#                               wants feedback on how many feature sets were removed.
#
# RETURNS
# 1. filtered_featureSet_array: The inputted feature-set array but with the prevlance
#                               violating feature sets removed.
# 2. fs_removed_lower:          A numpy array of the names of feature sets removed 
#                               because their prevalence was too low.
# 3. n_fs_removed_lower:        The count of feature sets removed for high low
#                               prevalence.
# 4. fs_removed_upper:          A numpy array of the namesof feature sets removed 
#                               because their prevalence was too high.
# 5. n_fs_removed_upper:        The count of feature sets removed for having high
#                               prevalence.
#

def boundaryfilter(my_featureSet_array, verbose = True):
    # Identify the feature sets that are too few, and extract the column names.
    fs_removed_lower = \
        my_featureSet_array.loc[:, 
        numpy.insert(
        ((my_featureSet_array.loc[:, my_featureSet_array.columns != 'person_id'] != 0).sum(axis=0) < min_criterion_count).values
            ,0
            ,False)
                   ].columns
    # Extract the count of feature sets that are too few.
    n_fs_removed_lower = len(fs_removed_lower)

    # Identify the feature sets that are too many, and extract the column names.
    fs_removed_upper = \
        my_featureSet_array.loc[:, 
        numpy.insert(
        ((my_featureSet_array.loc[:, my_featureSet_array.columns != 'person_id'] !=0).sum(axis=0) > max_criterion_count).values
            ,0
            ,False)
                   ].columns
    # Extract the count of feature sets that are too many.
    n_fs_removed_upper = len(fs_removed_upper)

    # Remove the feature sets that are no within the prevalence bounds.
    filtered_featureSet_array = \
        pandas.DataFrame(my_featureSet_array.drop(columns = numpy.insert(fs_removed_lower, 0, fs_removed_upper)))
    
    # Print message if arg{verbose} = True.
    if verbose == True:
        print("\n Filtering complete...")
        print("\t", len(filtered_featureSet_array.columns)-1, " feature sets remain.")
        print("\t", n_fs_removed_lower+n_fs_removed_upper, " feature sets removed, in total.")
        print("\t", n_fs_removed_lower, " feature sets removed because of low prevalence.")
        print("\t", n_fs_removed_upper, " feature sets removed because of high prevalence.")
    
    # Return outputs
    return [filtered_featureSet_array, fs_removed_lower, n_fs_removed_lower, fs_removed_upper, n_fs_removed_upper]


# -

# ### `evaloutputs()`: Compute the evaluation outputs

# + jupyter={"source_hidden": true} tags=[]
# A function to compute the evalaution outputs. The function automatically saves the
# contingency table and also returns it.
#
# ARGUMENTS
# 1. vec_featureSet:      A column from a pandas.Dataframe containing the feature set
#                         that needs evaluating.
# 2. vec_caseness:        A column from a pandas.Dataframe containing the caseness
#                         variable of interest.
# 3. savelocation:        The folder location where the output should be saved.
#
# RETURNS
# 1. prevalence:          The proportion of patients satisfying the definition of the
#                         feature set.
# 2. cba:                 Class balanced accuracy - the lower bound of the average
#                         sensitivity and average positive predictive value
#                         (a.k.a. precision) for all caseness values.
# 3. oddsRatio:           The ratio of the odds of caseness given the presence of feature
#                         set, to the odds of CMHD given the absence of the feature set.
#                         It can also be thought of as the multiplicative difference
#                         between correct and incorrect classification.
# 4. ppv:                 The proportion of patients satisfying the definition of the
#                         feature set who satisfy the caseness.
# 5. npv:                 The proportion of patients who do not satisfy the definition
#                         of the feature set who do not satisfy the caseness.
# 6. tn:                  The count of true negatives, i.e. the count of patients whose 
#                         feature-set value and caseness value are both zero.
# 7. fn:                  The count of false negatives, i.e. the count of patients whose 
#                         feature-set value is zero but whose caseness value is one.
# 8. fp:                  The count of false positives, i.e. the count of patients whose 
#                         feature-set value is one but whose caseness value is zero.
# 9. tp:                  The count of true positives, i.e. the count of patients whose 
#                         feature-set value and caseness value are both one.
#
def evaloutputs(vec_featureSet,
                vec_caseness):
    # ## Assess argument validty.
    
    # Check that both vectors are the same length.
    if len(vec_featureSet) != len(vec_caseness):
        print("\n**",
              "Feature-set and caseness vectors are of different lengths.",
             "**\n")
        return
    
    # ## Contingency table.
    # Make contingency table.
    contingencyTable = \
        pandas.crosstab(
            index = vec_featureSet,
            columns = vec_caseness
    )
    
    # Extract components of contingency table
    tn = contingencyTable.loc[0,0]
    fn = contingencyTable.loc[0,1]
    fp = contingencyTable.loc[1,0]
    tp = contingencyTable.loc[1,1]
    
    # ## Compute outputs.
    
    # Prevalence value per 1,000.
    #
    # I use 1 minus the prevalence of zeros because that
    # combines all the possibly-many values that indicate
    # the presence of the feature set.
    prevalence = \
        (1 - (sum(vec_featureSet == 0) / len(vec_featureSet))) * 10
    if prevalence < 0.01:
         prevalence = '< 0.01'
    else:
         prevalence = round(prevalence, 2)
    
    # Class balance accuracy.
    cba = \
        round( 0.5 * \
              ( (tp / max( (tp + fn), (tp + fp) ) ) + \
               (tn / max( (tn + fp), (tn +fn) ) ) ), 2)
    if cba < 0.01:
        cba = '< 0.01'
    
    # Odds ratio.
    if min( (tp * tn) , (fp * fn) ) == 0:
        oddsRatio = 'Undefined: One of the odds is zero.'
    else:
        oddsRatio = round( (tp * tn) / (fp * fn), 2)
    
    # Positive predictive value.
    ppv = 0.00 if (tp + fp) == 0 else tp / (tp + fp)
    if ppv > 0 and ppv < 0.01:
        ppv = '< 0.01'
    elif ppv < 1 and ppv > 0.999:
        ppv = '\u2248 1.00'
    else:
         ppv = round(ppv, 2)
         
    # Negative predictive value.
    npv = 0.00 if (tn + fn) == 0 else tn / (tn + fn)
    if npv > 0 and npv < 0.01:
        npv = '< 0.01'
    elif npv < 1 and npv > 0.999:
        npv = '\u2248 1.00'
    else:
         npv = round(npv, 2)
    
    
    
    return prevalence, cba, oddsRatio, ppv, npv, tn, fn, fp, tp


# -

# ### `evaleachcaseness()`: Iterate through the three caseness variables with a given feature set, and call the evalaution function.

# + jupyter={"source_hidden": true} tags=[]
# A function to iterate through the three caseness variables with a given feature set
# and call the evaluation function, evaloutputs.
#
# ARGUMENTS
# 1. vec_featureSet:      A column from a pandas.Dataframe containing the feature set
#                         that needs evaluating.
# 2. vec_caseness:        A column from a pandas.Dataframe containing the caseness
#                         variable of interest.
# 3. savelocation:        The folder location where the output should be saved.
#
# RETURNS
# n/a
#

def evaleachcaseness(vec_featureSet,
                     array_caseness,
                     savelocation = None):
    counter = 0
    for vec_caseness in array_caseness[["CMHD_dx_and_rx", "CMHD_rx_not_dx","CMHD_control"]]:
        contingencyTable,
        prevalence[counter],
        cba[counter],
        oddsRatio[counter],
        ppv[counter],
        npv[counter] = \
            evaloutputs(vec_featureSet,
                        vec_caseness,
                        savelocation = None)
        counter = counter + 1


# -

# ### `mutlinomRepresentation()`: Compute the multinomial representation of a feature set

# + jupyter={"source_hidden": true} tags=[]
# A function to compute the multinomial representation of a feature set.
#
# ARGUMENTS
# 1. var_vals:      An n-by-fs pandas.Dataframe of n patients and fs
#                   features, containing the feature sets that we want to
#                   compress into a single, multinomial representation.
#
# RETURNS
# 1. featureSet:    An n-by-1 pandas.Dataframe containing a multinomial
#                   representation of the inputted feature sets.
# 2. next_iter:     An indicator variable that is used by the parent
#                   function featuresetmi().
#
def mutlinomRepresentation(var_vals):
    # Check that the variables have more than three values and
    # only progress if False.
    for i_col in range(var_vals.shape[1]-1):
        unique_feature_vals = var_vals.iloc[:, i_col].drop_duplicates()
        if (len(unique_feature_vals) > 3):
            print("\n** Error: At least one of the",
                  "component features has more than",
                  "three values so the multinomial",
                  "representation will not be computed.**")
            print("Offending variable:", var_vals.columns.values[i_col],
                 ":", var_vals.columns.values[i_col].unique())
            unique_feature_vals
            next_iter = True
            return 0, next_iter

    # Get all combinations of values of the component features
    # and define feature set values for each multinomial combination.
    feature_combins = var_vals.drop_duplicates()
    feature_combins =\
        pandas.DataFrame(data = feature_combins, columns = var_vals.columns)\
        .reset_index()\
        .drop(['index'], axis = 1)
    feature_combins['multinom_vals'] = feature_combins.index


    # Define a vector indicating the feature set value.
    myMerge =\
        pandas.merge(
            var_vals,
            feature_combins,
            how = 'left',
            on = list(var_vals.columns.values)
    )

    # Extract multinomial representation as output variable.
    featureSet = myMerge['multinom_vals']
    next_iter = False
    return featureSet, next_iter


# + [markdown] tags=[]
# ### `featuresetmi()`: Calculate two-way mutual information for a features of order m.

# + jupyter={"source_hidden": true} tags=[]
# A function to calculate the two-way mutual information for a feature set.
#
# ARGUMENTS
# 1. featureSet_array:   An n-by-fs pandas.Dataframe of n patients and fs feature
#                        sets, or an fs-by-1 pandas.Dataframe containing the names
#                        of feature sets. If the fs-by-1 dataframe, then it is
#                        assumed that the feature sets are SNOMED-CT codes that 
#                        can be queried in the WHERE clause of BigQuery syntax.
# 2. casenessVector:     A column from a pandas.Dataframe containing the caseness
#                        variable of interest.
# 3. m:                  The order of feature set to be tested. 1 = Individuals, 
#                        2 = Pairs, 3 = Triplets.
# 4. savelocation:       The folder location where the output should be saved.
# 5. representation:     A choice of {'all', 'multi'} where 'all' = the feature set
#                        value is 1 when all components are 1, or 0 otherwise, and
#                        'multi' = the feature set values represent each combination
#                        of components' values.
# 6. source:             The source of the feature set: {'database', 'clinicial',
#                        'literature', 'interviews', 'PPI'}.
# 7. df_ppl_and_codes:   An optional argument that contains a list of all patients
#                        and all SNOMEDCT-CT codes that they have that are within
#                        the prevalence bounds.
# 8. verbose:            An optional argument to indicate whether the user wants
#                        feedback on how many feature sets were removed and saved.

#
# RETURNS
# n/a
#  
def featuresetmi(featureSet_array,
                 casenessVector,
                 m = None,
                 savelocation = None,
                 representation = None,
                 source = None,
                 df_ppl_and_codes = None,
                 verbose = False):
    # ## Assess argument validty.
    
    # Check order of feature set. If not provided,
    # default to m = 1.    
    if m == None:
        order_int = 1
        order_label = "Individuals"
        print("\nNo value for m provided." +
              "\n...Default value of m = 1 will be used.")
    elif m == 1:
        order_int = m
        order_label = "Individuals"
    elif m == 2:
        order_int = m
        order_label = "Pairs"
    elif m == 3:
        order_int = m
        order_label = "Triplets"
    else:
        print("\n** Error: Integer value between 1",
              "and 3 not supplied for m.**\n")
        return

    # Check and set save location.
    if savelocation == None:
        savelocation = \
           ("Mutual information saves/"+\
            order_label)
        print("\nNo save location provided." +
              "\n...Defaulting to ~/" + savelocation)    

    # ## Check encoding. If not provided, 
    # ## default to OR encoding.
    if representation == None:
        representation_label = "ALL"
        print("\nNo representation provided." +
              "\n...Defaulting to '" + representation_label + "' representation.")
    elif representation == "all":
        representation_label = "ALL"
    elif representation == "multi":
        representation_label = "MULTI"
    else:
        print("\n** Error: Representation value from ",
              "{'and', 'multi'} not provided.**\n")
        return
    
    # ## Check the source argument is provided.
    if source == None:
        print("\n** Error: No source argument provided.",
              "**\n")
        return
    
    # ## Set save string for particular caseness variable.
    caseness_type = casenessVector.columns.values[-1]
    if caseness_type == 'CMHD': 
        caseness_label = 'multinomial'
    elif caseness_type == 'CMHD_dx_and_rx': 
        caseness_label = 'definite'
    elif caseness_type == 'CMHD_rx_not_dx': 
        caseness_label = 'possible'
    elif caseness_type == 'CMHD_control': 
        caseness_label = 'control'
    
    print("\n\n\n****************************************")  
    print("Calculating mutual information values...")
    
    # Instantiate specific storage for mutual information.
    #global featureSet_MI
    featureSet_MI = \
        pandas.DataFrame(columns = ['Feature_set', 'Mutual_information'])

    # Instantiate batch number.
    #global batch 
    batch = 0

    # Instantiate tally of feature sets that are dropped due to low entropy.
    #global drop_tally
    drop_tally = 0

    # Define entropy of the particular caseness variable.
    #global entropy_caseness
    entropy_caseness = \
        scipy.stats.entropy(casenessVector.iloc[:,-1].value_counts(),
                            base = math.e)
    if verbose == True:
        print("Threshold f_MI is",entropy_caseness)
    
    # Check if the supplied feature set array is an n-by-fs array
    # of feature sets, or an fs-by-1 vector of feature-set names.
    # If it is the fs-by-1 vector, then pass all arguments to the
    # appropriate function; else, continue with the code.
    if featureSet_array.shape[1] == 1:
        # Display messages to user.
        if verbose == True:
            print("Overriding 'verbose == True' because there are too many feature sets.",
                 "Messages would become unweidly and slow down processing.")
            verbose = False
            # Inform users about limitation with database feature sets.
            if representation_label != "ALL":
                print("Only representation = 'all' is available for database feature sets.")
                return
        
        # The IPYNB file has already been run in this notebook but I'm repeating the run based
        # on guidance from this blog:
        # https://medium.com/@grvsinghal/speed-up-your-python-code-using-multiprocessing-on-windows-and-jupyter-or-ipython-2714b49d6fac
        # %run 'UNSEEN_helper_functions.ipynb'

        # Define a function to portion my iterable.
        #
        # https://stackoverflow.com/questions/51446327/python-3-generator-comprehension-to-generate-chunks-including-last
        def portion_maker(gen, portion_size):
            it = iter(gen)
            while True:
                portion = [*itertools.islice(it, 0, portion_size)]
                if portion:
                    yield portion
                else:
                    break
        
        # Do the main job of assessing the feature sets.
        if __name__ ==  '__main__': 
            gen = itertools.combinations(df_fs_database.snomedcode, m)
            portion_size = 7
            n_workers = 4
            for portion in portion_maker(gen, portion_size):
                print(f"This batch is {portion}.")
                list(Pool(n_workers).starmap(processdatabasefs, portion))
        
        # Calculate the total count of feature sets processed.
        count_fs = sum(1 for _ in itertools.combinations(featureSet_array.snomedcode, m))
    else:            
        # Ensure feature-set and casenesss values are matched for person_id.
        full_array = pandas.merge(casenessVector,
                                  featureSet_array,
                                  on = 'person_id',
                                  how = 'left')  

        # Define the m-way tuples of features sets as a numpy array. We will loop
        # through the rows of this array to create the feature sets.
        combins = \
            numpy.asarray(
                list(
                    itertools.combinations(
                        featureSet_array.columns[featureSet_array.columns != 'person_id'],
                        order_int)
                    )
                )

        # ## loop through the feature sets.
        for i_fs in tqdm(range(len(combins)), unit = 'feature set'):
            
            # Define an array indicating the feature set value.
            var_vals = full_array[combins[i_fs]]

            # Name the feature set.
            name_var = '-'.join(var_vals.columns.values)

            if representation_label == "ALL":
                fs_val = var_vals.all(True)
            elif representation_label == "MULTI":
                fs_val, next_iter = mutlinomRepresentation(var_vals)
                if next_iter:
                    drop_tally += 1
                    if verbose == True:
                        print("Dropped", name_var, "because at least one",
                              "component features has more than three values")
                    continue


            # Calculate the mutual information between the feature set and
            # caseness variable.
            f_MI = sklearn.metrics.mutual_info_score(fs_val, full_array[caseness_array.columns[1]])

            if f_MI < entropy_caseness:
                drop_tally += 1
                if verbose == True:
                    print("Dropped",name_var,"because f_MI =", f_MI)
                continue
            else:            
                # Store the name and mutual information value.
                if verbose == True:
                    print("\tSaved",name_var,"because f_MI =", f_MI)
                featureSet_MI.loc[len(featureSet_MI),:] = name_var, f_MI

            if len(featureSet_MI) > 9:
                    # Increment batch.
                    batch += 1

                    # Make an interim save of results.
                    current_dir = os.getcwd()
                    print(current_dir + savelocation)
                    os.chdir(current_dir + "/" + savelocation)
                    pyarrow.parquet.write_table(pyarrow.Table.from_pandas(featureSet_MI),   
                                                source + "_" +
                                                caseness_label + "_" +
                                                representation_label + "_" + 
                                                "batch" + 
                                                str(batch) + 
                                                ".parquet")
                    os.chdir(current_dir)
                    featureSet_MI.to_csv(savelocation + "/" +
                                         source + "_" +
                                         caseness_label + "_" + 
                                         representation_label + "_" + 
                                         "batch" + 
                                         str(batch) + 
                                         ".csv", index = False)
                    print("\nInterim save made")
                    # Instantiate new storage.
                    featureSet_MI = \
                        pandas.DataFrame(columns = ['Feature_set', 'Mutual_information'])
                    
                    # Calculate the total count of feature sets processed.
                    count_fs = len(combins)
                    
    # Increment counter.
    batch += 1
    # Final save.
    if len(featureSet_MI) != 0:
        current_dir = os.getcwd()
        print(current_dir + savelocation)
        os.chdir(current_dir + "/" + savelocation)
        pyarrow.parquet.write_table(pyarrow.Table.from_pandas(featureSet_MI),
                                    source + "_" +
                                    caseness_label + "_" +
                                    representation_label + "_" + 
                                    "batch" + 
                                    str(batch) + 
                                    ".parquet")
        os.chdir(current_dir)
        featureSet_MI.to_csv(savelocation + "/" +
                             source + "_" +
                             caseness_label + "_" +
                             representation_label + "_" +
                             "batch" + 
                             str(batch) + 
                             ".csv", index = False)

    # Feedback messages.
    print("...\n")
    print(str(batch), "batch(es) of feature sets processed.")
    print(str(drop_tally), "/",
          str(count_fs),
          "feature sets dropped due to low entropy.")
    print("****************************************")  


# -

# ### `processdatabasefs()`: Calculate two-way mutual information for *database* features sets, specifically.

# + jupyter={"source_hidden": true} tags=[]
def processdatabasefs(*snomedcodes):
    global df_ppl_and_codes
    global caseness_array
    global drop_tally
    global batch
    global featureSet_MI
    global entropy_caseness
    
    # Define the feature set.
    fs = [snomedcodes]
    
    # Name the feature set.
    name_var = '-'.join(map(str,list(fs)))

    # Check to see if anyone in the cohort has the codes in their record.
    df_temp = df_ppl_and_codes
    df_temp['ary_pAc_in_fsdf'] = df_ppl_and_codes.snomedcode.isin(fs)
    do_patients_have_all_query_codes = \
        (df_temp.groupby(['person_id']).ary_pAc_in_fsdf.count() == len(fs)).astype(int)

    # Join feature set with caseness_array to match person_id.
    fs_val = pandas.merge(caseness_array,
                          do_patients_have_all_query_codes,
                          how = 'left',
                          on = 'person_id')

    # Calculate the mutual information between the feature set and
    # caseness variable.
    f_MI = sklearn.metrics.mutual_info_score(fs_val['ary_pAc_in_fsdf'], fs_val[caseness_array.columns[1]])

    if f_MI < entropy_caseness:
        print(f'f_MI < entropy_caseness')
        drop_tally += 1
    else:            
        # Store the name and mutual information value.
        featureSet_MI.loc[len(featureSet_MI),:] = name_var, f_MI

    if len(featureSet_MI) > 9:
            # Increment batch.
            batch += 1

            # Make an interim save of results.
            current_dir = os.getcwd()
            print(current_dir + savelocation)
            os.chdir(current_dir + "/" + savelocation)
            pyarrow.parquet.write_table(pyarrow.Table.from_pandas(featureSet_MI),   
                                        source + "_" +
                                        caseness_label + "_" +
                                        representation_label + "_" + 
                                        "batch" + 
                                        str(batch) + 
                                        ".parquet")
            os.chdir(current_dir)
            featureSet_MI.to_csv(savelocation + "/" +
                                 source + "_" +
                                 caseness_label + "_" + 
                                 representation_label + "_" + 
                                 "batch" + 
                                 str(batch) + 
                                 ".csv", index = False)
            print("\nInterim save made")
            # Instantiate new storage.
            featureSet_MI = \
                pandas.DataFrame(columns = ['Feature_set', 'Mutual_information'])


# -

# ### `portion_maker()`: Define a function to portion my iterable.
# https://stackoverflow.com/questions/51446327/python-3-generator-comprehension-to-generate-chunks-including-last
def portion_maker(gen, portion_size):
    it = iter(gen)
    while True:
        portion = [*itertools.islice(it, 0, portion_size)]
        if portion:
            yield portion
        else:
            break


# ### `init_worker()`: Initialize worker processes for parallell processing.
def init_worker(par_ppl_and_codes, par_caseness_array, par_drop_tally, par_batch, par_featureSet_MI, par_entropy_caseness):
    global df_ppl_and_codes
    global caseness_array
    global drop_tally
    global batch
    global featureSet_MI
    global entropy_caseness
    # store argument in the global variable for this process
    df_ppl_and_codes = par_ppl_and_codes
    caseness_array = par_caseness_array
    drop_tally = par_drop_tally
    batch = par_batch
    featureSet_MI = par_featureSet_MI
    entropy_caseness = par_entropy_caseness
    #print(f'df_ppl_and_codes.head() = {df_ppl_and_codes.head()}')


# ### patienthassnomed(): Check whether a patient's record includes some particular SNOMED-CT codes.

# + jupyter={"source_hidden": true} tags=[]
def patienthassnomed(fs_df, snomedcode):
    ary_pAc_in_fsdf = []
    if math.isnan(i):
        ary_pAc_in_fsdf.append(False)
    else:
        ary_pAc_in_fsdf.append((i == fs_df).any()[0])
    return ary_pAc_in_fsdf


# + [markdown] tags=[]
# ### featuresetmi_database(): Calculate two-way mutual information for a *database* features of order m.
# -

# ### getfsarray(): Produce a database feature set for submission to `featuresetmi()` function

# + jupyter={"source_hidden": true} tags=[]
def getfsarray(fs_df):
    
    # Instantiate a Google BigQuery client.
    client = bigquery.Client()
    
    # Define the BigQuery syntax.

    sql_CTEs_body = \
    """
    # Make a table of person ID and their SNOMED-CT codes from the list of codes of interest.
    WITH tbl_persons_and_codes AS
    (
    SELECT
        DISTINCT person_id
        ,src_snomedcode
    FROM
        yhcr-prd-phm-bia-core.CY_FDM_PrimaryCare_v5.tbl_SRCode
    WHERE
        src_snomedcode IN ('""" + '\', \''.join(map(str, fs_df['src_snomedcode'].to_list())) + """')
    )
    """

    sql_pivot = \
    """
        SELECT
        CONCAT("SELECT person_id,"
               , STRING_AGG(CONCAT("CASE WHEN src_snomedcode='",src_snomedcode,"' THEN 1 ELSE 0 END AS `_",src_snomedcode,"`")), 
            " FROM `tbl_persons_and_codes`",
            " GROUP BY person_id, src_snomedcode ORDER BY person_id")
    FROM (  SELECT DISTINCT src_snomedcode FROM `tbl_persons_and_codes` ORDER BY src_snomedcode  )
    """
    sql = client.query(sql_CTEs_body + sql_pivot).to_dataframe()['f0_'].iloc[0]
    
    feastureSet_array = \
        client.query(sql_CTEs_body +
                     sql).to_dataframe()
    
    return feastureSet_array


# + [markdown] tags=[]
# ### NOT IN USE calculatemi(): Do the work of calculating the MI between a given feature set and caseness variable

# + jupyter={"source_hidden": true} tags=[]
def calculatemi(var_vals,
                name_var,
                representation_label,
                drop_tally,
                verbose):
    
    # Formulate the representation.
    if representation_label == "ALLrepresentation":
        fs_val = var_vals.all(True)
    elif representation_label == "MULTIrepresentation":
        fs_val, next_iter = mutlinomRepresentation(var_vals)
        if next_iter:
            return next_iter

    # Calculate the mutual information between the feature set and
    # caseness variable.
    f_MI = sklearn.metrics.mutual_info_score(fs_val, full_array.iloc[:,-1])

    if f_MI < entropy_caseness:
        drop_tally += 1
        if verbose == True:
            print("Dropped",name_var,"because f_MI =", f_MI)
    else:            
        # Store the name and mutual information value.
        if verbose == True:
            print("\tSaved",name_var,"because f_MI =", f_MI)
        featureSet_MI.loc[len(featureSet_MI),:] = name_var, f_MI

    if len(featureSet_MI) > 9:
            # Increment batch.
            batch += 1

            # Make an interim save of results.
            featureSet_MI.to_csv(savelocation +
                              representation_label + "_" +
                              "_batch" + \
                              str(batch) + \
                              ".csv", index = False)
            print("\nInterim save made")
            # Instantiate new storage.
            featureSet_MI = \
                pandas.DataFrame(columns = ['Feature set', 'Mutual information'])
