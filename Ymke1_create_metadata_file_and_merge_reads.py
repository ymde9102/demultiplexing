
# coding: utf-8

# In[3]:

import pandas

primer_sheet = pandas.read_csv("Primers_and_pairs_plus_switched.csv", header=5, skip_blank_lines=5)
list(primer_sheet)


# Need to drop some of the columns as they are specifications for other experiments.

# In[28]:

primer_sheet = primer_sheet.drop({'Unnamed: 3','Unnamed: 7', 
                                  'Unnamed: 12', 'Unnamed: 14', 
                                  'Forward primer.1', 'Reverse primer.1',
                                  'Sample.1'}, axis=1)
primer_sheet.rename(columns={'Forward primer': 'f_primer_ymke1',
                            'Reverse primer' : 'r_primer_ymke1',
                            'Sample' : 'sample_ymke1'}, inplace=True)


# In[29]:

list(primer_sheet)


# In[30]:

primer_dict = dict(zip(primer_sheet['Primer name'], primer_sheet['Complete sequence.']))
primer_dict.viewkeys()


# ### Create Matrix file as required by QIIME

# There is a page at [QIIME](http://qiime.org/1.7.0/documentation/file_formats.html) website that helps deal with this. It is not tolerant of deviations, so one has to make sure it is properly formatted ( it can only check so many things though ).

# The forward and reverse primer information then needs to be mapped to the keys juset generated. 
# 
# If one looks at the original sheet, there is a description of the primers used and their sequences, as well as how these were then used in the sequencing experiments. Ideally, these should have been two files, a reference file to show what was used, and another to show how it was used. It is often a good idea to separate data like this.

# In[31]:

primer_sheet['ReversePrimer'] = primer_sheet['r_primer_ymke1'].map(primer_dict)


# In[32]:

primer_sheet['LinkerPrimerSequence'] = primer_sheet['f_primer_ymke1'].map(primer_dict)


# In[33]:

print(primer_sheet[0:5])


# In[34]:

primer_sheet = primer_sheet.drop({'Primer', 'Primer name', 'Complete sequence.', 'Barcode', 'Spacer'}, axis=1)
# drop the empty rows
primer_sheet = primer_sheet.dropna(how='all')


# In[35]:

# see how big it is
print(primer_sheet[48:55])


# Use the `LinkerPrimerSequence` to extract the barcode + spacer sequence for each sample.

# In[36]:

for index, row in primer_sheet.iterrows():
    new_value = row['LinkerPrimerSequence'][0:7]
    primer_sheet.loc[index, 'BarcodeSequence'] = new_value    


# In[37]:

primer_sheet[45:55]


# In[38]:

for index, row in primer_sheet.iterrows():
    new_value = row['ReversePrimer'][0:7]
    primer_sheet.loc[index, 'ReverseBarcodeSequence'] = new_value


# In[39]:

primer_sheet = primer_sheet.rename(columns={'sample_ymke1': '#SampleID'})
print(primer_sheet[45:55])


# QIIME is particular about what characters are in the SampleID column.

# In[40]:

list(primer_sheet)


# In[41]:

# reorder the columns so the mapping checker doesn't complain about formatting
new_index = ['#SampleID',
 'BarcodeSequence',
 'LinkerPrimerSequence',
 'f_primer_ymke1',
 'r_primer_ymke1', 
 'ReversePrimer', 
 'ReverseBarcodeSequence']

primer_sheet = primer_sheet.reindex(columns=new_index)


# In[42]:

primer_sheet['#SampleID'] = primer_sheet['#SampleID'].str.replace(':', '.')
primer_sheet['Description'] = primer_sheet['#SampleID']


# In[43]:

print(primer_sheet[0:5])
list(primer_sheet)


# ### Change the names of duplicated samples

# In[44]:

repeated_records_list = ['1.12','2.12','3.12','4.12']

records_2b_changed_dict = {'1.13' : ['IllB_48_L_341F_BA','IllA_54_L_805R_BA'], 
                      '2.13' : ['IllB_45_L_341F_BA','IllA_53_L_805R_BA'],
                      '3.13' : ['IllB_43_L_341F_BA','IllA_52_L_805R_BA'],
                      '4.13' : ['IllB_42_L_341F_BA','IllA_51_L_805R_BA']}
records_2b_changed_dict.viewvalues()


# In[45]:

for index, row in primer_sheet.iterrows():
        if row['#SampleID'] in repeated_records_list:
            #print(row['#SampleID'])
            for key, values in records_2b_changed_dict.iteritems():
                if values[0] == row['f_primer_ymke1'] and values[1] == row['r_primer_ymke1']:
                    #print(values[0], row['f_primer_ymke1'], values[1], row['r_primer_ymke1'])
                    primer_sheet.set_value(index, '#SampleID', key)


# In[46]:

primer_sheet['Description'] = primer_sheet['#SampleID']


# In[47]:

mapping_filename = "barcode_sequences_ymke1.tsv"
primer_sheet.to_csv(mapping_filename, sep="\t", index=False)


# ### Validate the metadata mapping file

# Validate the mapping file to make sure it is compatible with QIIME with the QIIME script.

# In[48]:

get_ipython().system(u'validate_mapping_file.py -m barcode_sequences_ymke1.tsv')


# Errors and warnings pertain to duplicate barcode sequences. Of course it is not a duplicate when the primer sequence is included, and hopefully the next step handles this. Formatting seems fine otherwise.

# ### Merging the FASTQ files

# We can try the program [PEAR](http://sco.h-its.org/exelixis/web/software/pear/). Changed the path to reflect it on your system.

# In[22]:

get_ipython().system(u'/pear-0.9.10-bin-64/pear-0.9.10-bin-64 -j 2 -f Ymke1_S1_L001_R1_001.fastq -r Ymke1_S1_L001_R2_001.fastq -o Ymke1_S1_L001_merged.fastq')


# Seems like this one had no issues...perhaps the QIIME script is not as robust. This merge takes a long time to run which is not reflected in the output.
