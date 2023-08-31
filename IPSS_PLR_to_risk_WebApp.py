##########################################################################################
##### Dec 2022 Roche Hackathon                                                    ########
##### Team: EpiMeter                                                              ########
##### Author: Julie Huang                                                         ########
##### Description: Predict effect of oligo mismatches on assays                   ########
##### As an application of Intelligent pathogen surveillance system, this tool    ########
##### leverages mismatch information from PLR (PCR assay post launch report), to  ########
##### predict emerging isolate sequences at high risk of being missed by the      ########
##### current launched assay design. The results can be launched as a WebApp by   ########
##### Streamlit tool                                                              ########            
##########################################################################################

import os
import pandas as pd
import numpy as np
import streamlit as st

oligo_nt_set = {'a', 'c', 'g', 't', '-', 'A', 'C', 'G', 'T'}

def count_no_mismatches(oligo_string):
    '''
    Calculate no. mismatch and frac of mismatch in the oligo
    '''
    no_match_count = 0
    for letter in oligo_string:
        if letter != "." and letter !="-":
            no_match_count+=1
    fraction_mismatch = round(no_match_count/len(oligo_string), 2)
    return (no_match_count, fraction_mismatch)


def three_prime_end_mismatch(oligo_string, no_3_prime_nts):
    '''
    Calculate no. of mismatch at 3' end of the oligo
    '''
    non_match_count = 0
    three_prime = oligo_string[-no_3_prime_nts:]
    for letter in three_prime:
        if letter != ".":
            non_match_count+=1
    return non_match_count


def count_primer_relative_freq(oligo_sequence, df, oligo_column_name):
    df_count_agg = df[['id', oligo_column_name]].groupby(oligo_column_name).count().reset_index()
    freq_count = df_count_agg.loc[df_count_agg[oligo_column_name]==oligo_sequence, 'id'].values[0]
    total_ids = df.shape[0]
    freq_as_ratio = round(int(freq_count)/int(total_ids), 2)
    return freq_count, freq_as_ratio


def input_PLR_return_taxonIDs_need_notice(PLR_file, PLR_sheet_name, no_top_taxon_to_show, three_prime_distance):
    '''
    PLR_file: path to the PLR excel sheet
    PLR_sheet_name: the name of the excel tab of "*_report", with oligo related information
    no_top_taxon_to_show: no. of the most worth noticing taxons to show, that may not be amplified by assay
    three_prime_distance: for forward and reverse primers, the n nucleotides to 3' end that needs mismatch check
    '''
    ## Read in dataframe, oligo names and column names
    df = pd.read_excel(PLR_file, sheet_name=PLR_sheet_name)
    oligo_header = df.columns.to_list()
    complete_header = df.iloc[0].to_list()

    ## Merge and organize new header
    new_header = []
    for column_no in range(len(oligo_header)):
        if oligo_header[column_no].startswith("Unnamed"):
            new_col_name = complete_header[column_no]
        elif set(complete_header[column_no]).issubset(oligo_nt_set):
            new_col_name = oligo_header[column_no]
        else:
            new_col_name = str(str(oligo_header[column_no]).split(".")[0]) + ":" + str(complete_header[column_no])
        new_header.append(new_col_name)

    ## Rename df with new_header
    df = df[1:]
    df.columns = new_header

    ## For primer oligo sequence columns ending with (f) and (r), make mismatch related calculations
    ## For primer oligo sequence, also count frequency and relative abundance of the haplotype
    primer_sequence_columns = []
    three_prime_rank_columns = []
    no_mismatch_rank_columns = []
    mismatch_ratio_rank_columns = []
    haplotype_count_columns = []
    for col_names in new_header:
        #print("col_names", col_names)
        if col_names.endswith("(f)") or col_names.endswith("(r)"): # the oligo columns
            three_prime_mismatch_name = str(col_names) + "_3_prime_mismatch"
            mismatch_name = str(col_names) + "_mismatch_count"
            mismatch_ratio_name = str(col_names) + "_mismatch_ratio"
            haplotype_count_name = str(col_names) + "_haplotype_count"
            haplotype_ratio_name = str(col_names) + "_relative_abundance"

            df[three_prime_mismatch_name] = df[col_names].apply(three_prime_end_mismatch, args=(three_prime_distance, ))
            df[[mismatch_name, mismatch_ratio_name]]= df[col_names].apply(count_no_mismatches).tolist()
            df[[haplotype_count_name, haplotype_ratio_name]] = df[col_names].apply(count_primer_relative_freq, args=(df, col_names)).tolist()

            three_prime_rank_columns.append(three_prime_mismatch_name)
            no_mismatch_rank_columns.append(mismatch_name)
            mismatch_ratio_rank_columns.append(mismatch_ratio_name)
            haplotype_count_columns.append(haplotype_ratio_name)
            primer_sequence_columns.append(col_names)
            
            columns_for_strain_ranking = three_prime_rank_columns + \
            no_mismatch_rank_columns + mismatch_ratio_rank_columns + haplotype_count_columns
        ## Pending 
        ## need do develop for "p", once quencher and dye location obtained

    ascending_list = [False] * len(columns_for_strain_ranking)
    df_need_notice = df.sort_values(columns_for_strain_ranking, ascending = ascending_list)

    ## Rank strain IDs
    to_show_columns = ['Isolate_Id', 'Isolate_Name', 'Location', 'Host', 'Collection_Date', 'country']
    to_show_extended = to_show_columns + primer_sequence_columns + haplotype_count_columns \
    + three_prime_rank_columns + mismatch_ratio_rank_columns
    df_need_notice = df.sort_values(columns_for_strain_ranking, \
                                    ascending = ascending_list).head(no_top_taxon_to_show)[to_show_extended]  
    return df_need_notice



############################################################
## Streamlit Webapp implementation on an example dataset 
## Input an excel sheet containing PLR results

PLR_file = '/Users/huangz36/Documents/Hackathon_2022/flu_hackathon_1100.xlsx'
sheet_name = 'FluA_report'


#Back-up demo input data with more sequences curated
#PLR_file = '/Users/huangz36/Documents/Hackathon_2022/flu_hackathon_10000.xlsx'
#sheet_name = 'flu_hackathon_10000'



## Streamlit WebApp flow
st.title('IPSS: predict effects of oligo mismatches on assay performance')

@st.cache
def load_PLR(PLR_file, sheet_name):
	PLR_df = pd.read_excel(PLR_file, sheet_name=sheet_name)
	return PLR_df

data_load_state = st.text('Loading PLR...')
PLR_raw = load_PLR(PLR_file, sheet_name)
data_load_state.text('Loading PLR...done!')

st.subheader('Peek of raw PLR')
st.write(str(PLR_raw.shape[0]), " unique Seq IDs were included in the PLR")
st.write(PLR_raw.head(20))

st.subheader('Emerging Seq ID containing oligo mismatches may have the risk to be missed by current assay')
st.text('Enter parameters for risky sequence prediction')
taxon_no = st.text_input("How many sequences at risk do you want to see?", key="no_top_taxon_to_show")
three_prime = st.text_input("For designed primer, how close a mismatch to the 3' end (bp) do you care most?\n(Value must be >=1)", key="three_prime_distance")


if taxon_no and three_prime:
	st.text('Predicting Seq IDs at risk of being missed by assay...')
	output_risk_df = input_PLR_return_taxonIDs_need_notice(PLR_file, sheet_name, int(taxon_no), int(three_prime))
	st.table(output_risk_df)
	st.text('Please inspect the date, prevalence and country origin of the Seq ID as risk. \nTake action when necessary...')


output_for_slides = input_PLR_return_taxonIDs_need_notice(PLR_file, sheet_name, 100, 4)
output_for_slides.to_csv("/Users/huangz36/Documents/Hackathon_2022/flu_1100_output.txt", sep='\t', header=True, index=False)