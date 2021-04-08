## Written by: Kate Johnson, kej310@nyu.edu
"""
Adding amino acid and gene information

Requires the UPDATED snplist files and csv file with features interested in

Snplist files are generated with the readreport.py script.
Updated snplist files are generated with the ConsensusFasta.Coverage.v4.py script.
Amino acid snplist files are generated with the AddAminoGene.5.py script.

The csv file with features that are of interest should have the following columns:
SEGMENT,START,END,NAME

python3 AddAminoGene.5.py -- freqcut 0.001
                           --ref ../MergedRuns/reference/SARS-COV-2.fasta
                           --var ../MergedRuns/NYU/fullvarlist/
                           --strain COV19
                           --features ./sars-cov-2-features.3.csv
                           --save_dir ../MergedRuns/NYU/aasnplist/

"""
import os
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r', required=True,help='Indicate path and reference fasta file') #args.ref
parser.add_argument('--var','-v', required=True,help='Indicate path to variant files')
parser.add_argument('--strain','-s', required=True,help='Indicate strain')
parser.add_argument('--features','-f', required=True,help='Indicate csv file with features of interest')
parser.add_argument('--save_dir','-sd', required=True,help='Indicate where new file should save')
parser.add_argument('--freqcut','-q', default=0.001,help='Indicate frequency cutoff')
args = parser.parse_args()

#amino acid dictionary
aminoacid = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}

def list_duplicates(ntpos_list):
    """
    INPUT: Tims code will sometimes duplicate positions that have only 1 read
    input list of ntpos to see if there are duplicates
    OUTPUT: return list of positions that have been duplicated
    """
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set( x for x in ntpos_list if x in seen or seen_add(x) )
    print(list(seen_twice))
    return list( seen_twice )

def read_fasta(fp):
    """
    INPUT: Fasta file reference
    """
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def open_fasta(filename):
    """
    INPUT: fasta file reference
    OUTPUT: dictionary with segment name as key and sequence of segment as value
    """
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict

def creatList(r1,r2):
    """
    INPUT: the start and end of a gene position
    OUTPUT: list of positions within the range of the start and end of a gene
    """
    return list(range(r1,r2+1)) #python up to number you don't want included

def getindex(codon):
    """
    Identifying where the nucleotide falls in codon
    INPUT: major codon from snplist file
    OUTPUT: Index location of nucleotide
    """
    return [i for i, c in enumerate(codon) if c.isupper()]

def AdjustMinor(minordf):
    """
    INPUT: a codon df generated from translate function
    OUTPUT: a minor codon, minor amino acid, and nonysn info if present
    """
    print("Gathering minor information")
    minordf = minordf.reset_index(drop=True)
    for index, row in minordf.iterrows():
        codon = row['majorcodon']
        aa = row['majoraa']
        #print(codon)

        if pd.isna(codon):
            continue

        else:
            if '-' in codon: #if there is a deletion
                continue
            elif 'n' in codon or 'N' in codon: #if there is an N- will happen with primer trimming
                continue

            else:
                upperIndex = getindex(codon)[0] #identify where nucleotide is
                #print(upperIndex)
                minorcodon = list(codon) #make codon a list of three nucleotides

                minorcodon[upperIndex] = row['minor'].upper() #put minor nt in index position and make uppercase
                minorcodon = ''.join(minorcodon) #join into string (no longer list)
                #print(codon, minorcodon)

                if '-' in minorcodon: #if there is a deletion for minor nt
                    continue

                else: #grab aa info for min
                    minoraa = aminoacid[minorcodon.lower()] #amino acid letter
                    if minoraa == aa:#if minor aa is the same as the major amino acid
                        nonsyn = 'syn'
                    elif minoraa != aa:
                        nonsyn = 'nonsyn'

                    #add all info to the main dataframe
                    minordf.loc[index,"minorcodon"] = minorcodon
                    minordf.loc[index,"minoraa"] = minoraa
                    minordf.loc[index,"nonsyn"] = nonsyn

                    #print(minordf.loc[index])

    return minordf # return main dataframe of minor information



def translate(region_df,regionList,gene_name):
    """
    INPUT: subseted dataframe for region that we are interested
    OUTPUT: outputs a dataframe with updated information such as the major
    codon, amino acid, and gene id used with the region
    """
    print("Translating codons")

    region_df["major"] = region_df["major"].fillna("N") # change any missing data to N
    region_df = region_df.sort_values(by=["ntpos"])
    regionLength = len(regionList) #length of gene - must be divisible by three to translate
    aalist = list(range(0, regionLength, 3))  #generate a list going multiples of three
    dfntposLen = len(list(region_df['ntpos'])) #incase pos repeated due to poor alignments
    ntpos_list = list(region_df['ntpos'])

    #if regionLength % 3 == 0: #if no remainder after dividing by 3
    if dfntposLen % 3 == 0:
        print(dfntposLen%3)
        full_df = [] #empty df to append to
        for i in range(0, regionLength, 3): #iterate by 3's
            aaposition = aalist.index(i) + 1 #amino acid position - add 1 as python is 0based
            cod_df = region_df[i:i+3] #grab rows that correspond to amino acid
            cod_df = cod_df.reset_index(drop=True) #reset index (now 0-2)
            cList = list(region_df[i:i + 3]['major']) #grab the major nts for codon
            #print(cList)

            codonList = [x.lower() for x in cList] #make it lower to grab amino acid in dict
            codon = (''.join(codonList)).lower() #join the list of major nts to make string

            rlist = list(region_df[i:i + 3]['refnt'])
            refcodlist = [x.lower() for x in rlist]
            rcodon = (''.join(refcodlist)).lower()

            cod_df['gene_id'] = gene_name #add gene name for subset data
            cod_df['aapos'] = aaposition

            if 'n' in codon or '-' in codon:
                cod_df['gene_id'] = gene_name #add gene name for subset data
                cod_df['aapos'] = aaposition #add position
                full_df.append(cod_df) #but don't translate


            elif len(codon)==3:

                refaa = aminoacid[rcodon]
                cod_df['refcodon'] = rcodon
                cod_df['refaa'] = refaa

                AA = aminoacid[codon] #grab the amino acid letter using the codon string
                c = codonList
                c[0] = c[0].upper() #for first nt position in codon make it upppercase
                first = "".join(c) #join as string
                #print(type(first))
                #print(first)

                c[0] = c[0].lower() #change the first nt to lowercase
                c[1] = c[1].upper() #make second position in codon uppercase
                second = "".join(c)

                c[1] = c[1].lower() #make second position lower case
                c[2] = c[2].upper() #make third position uppercase
                third = "".join(c)

                cod_df['majoraa'] = AA #add major amino acid information to dataframe
                cod_df.loc[0,"majorcodon"] = first
                cod_df.loc[1,"majorcodon"] = second
                cod_df.loc[2,"majorcodon"] = third

                #cod_df['gene_id'] = gene_name #add gene name for subset data
                #cod_df['aapos'] = aaposition
                #cod_df.at[0,'majorcodon'] = str(first) #add codon for first row of subset data
                #cod_df.at[1,'majorcodon'] = second #add codon for second row of subset data
                #cod_df.at[2,'majorcodon'] = third #add codon for third row of subset data

                #print(cod_df)

                full_df.append(cod_df) #append to empty list

        full_df=pd.concat(full_df,sort=False) #concatenate empty list into full dataframe
        return full_df


    #elif regionLength % 3 != 0: #if not divisible by 3, still want to add info
    elif dfntposLen % 3 != 0: #if not divisible by 3, still want to add info
        print("Not divisible by 3 {0}".format(regionLength % 3))
        print(gene_name)
        x = list_duplicates(ntpos_list)
        cod_df = region_df #just add gene id column name
        cod_df = cod_df.assign(gene_id = gene_name)
        return cod_df

def NoGen(subsetDF):
    """
    if there isn't a gene id/name in a region (5'utrs/3'utrs) just add empty gene id
    to the dataframe
    """
    subsetDF= subsetDF.assign(gene_id="")
    return subsetDF

def markRegions(feature_df,df):
    """
    INPUT: Input a features dataframe and snplist variant file
    OUTPUT: Output info to add to new snplist file with aa information
    """
    print("Subsetting dataframe by specified region")
    total_df = [] #empty df list to append and concatenate later
    tot_Regions = [] #looking at all ntpos covered through iterations
    for index, row in feature_df.iterrows(): #iterate through featers df
        region = creatList(row['START'],row['END']) #generate list of ntpos
        tot_Regions = tot_Regions + region #add ntpos to list
        df_region = df[df.ntpos.isin(region)] #subset dataframe for ntpos only in region of interest
        #print(df_region)
        #check to make sure the df_region is the same length as the region we are interested in
        print(len(region))
        print(len(list(df_region['ntpos'])))

        gene_id = row['NAME'] #name of gene
        print(gene_id)
        if len(region) == len(list(df_region['ntpos'])):
            trans_df = translate(df_region,region,gene_id) #pipe through translate function to add info
            total_df.append(trans_df) #append the output df to total df

    tot_Regions = list(set(tot_Regions)) #take the set of the ntpos list to remove duplicates
    noRegion = NoGen(df[~df.ntpos.isin(tot_Regions)]) #tilde indicates not in

    total_df.append(noRegion)
    total_df=pd.concat(total_df, sort=False)
    total_df.sort_values(by=['ntpos'])
    return total_df

def AddRefSeq(df,refseqdf):
    """
    INPUT: Reference sequence
    OUTPUT: Data frame with reference nt attached
    """
    print("Adding Reference Seq")

    m = pd.merge(df, refseqdf, on="ntpos", how="outer")
    #m = m.reset_index() #reset index (now 0-2)

    return m
    #pd.merge(restaurant_ids_dataframe, restaurant_review_frame, on='business_id', how='outer')


    """
    for index, row in df.iterrows():
        ntpos_i = row['ntpos'] - 1
        refnt = refseq[ntpos_i]
        df.loc[index,"refnt"] = refnt

    return df
    """



### NEED TO GENERATE A CHECK THAT MAKES SURE ALL NTPOS IN THE REF ARE PRESENT AND IF NOT ADD WITH AN 'N'

#### RUN CODE BELOW ####
path = args.var
feat = pd.read_csv(args.features) #pandas will read the features csv file
STRAIN = args.strain
ref_dict = open_fasta(args.ref)
freqcut = args.freqcut
NT_LIST = ['A','T','G','C']



for SEGMENT in ref_dict:
    print(SEGMENT)
    seg_feats = feat[feat.SEGMENT == SEGMENT] #filter features for segment
    #print(seg_feats)
    refseq = list(ref_dict[SEGMENT]) #make the ref seq a list
    SegLength = list(range(1, len(refseq) + 1)) #generate list of ntpos for ref
    refd = {"refnt":refseq,"ntpos":SegLength} #build dict to build df
    refdf = pd.DataFrame(refd) #build ref df with ntpos and nt
    #print(refdf)


    for infile in glob.glob( os.path.join(path, '*{0}.{1}.Updated.{2}.snplist.csv'.format(STRAIN.upper(),SEGMENT,freqcut)) ):
        print(infile)
        fname = str(infile)
        #print(fname)
        df = pd.read_csv(infile) #pandas will read the csv file create dataframe
        df = df[df.segment == SEGMENT]
        #df = AddRefSeq(df,refseq)

        #drop duplicate ntpos from snplist files (happens if low coverage)
        df = df.drop_duplicates(subset=['ntpos'],keep='first')
        df = AddRefSeq(df,refdf) #add ref info to dataframe
        #name = list(set(list(df['name'])))#[0]
        #print(name)
        name = fname.split(path)[1].split(".")[0]
        print(name)


        UpdatedDF = markRegions(seg_feats,df)
        #print(UpdatedDF)

        dfmin = UpdatedDF[(UpdatedDF.minor.isin(NT_LIST))]
        nomin = UpdatedDF.merge(dfmin, how = 'outer' ,indicator=True).loc[lambda x : x['_merge']=='left_only']
        minchange = AdjustMinor(dfmin)
        #print(minchange)


        frames = [nomin,minchange]
        print('frames')
        #print(frames)

        print('----------------')
        print('concat')
        final_df = pd.concat(frames,sort=False)
        #print(final_df)

        print('----------')
        print('dropthings')
        #final_df=final_df.drop(columns=["index"])
        final_df = final_df.drop(columns=['_merge'])
        #final_df = final_df.drop(columns=['level_0']) #where does this come
        #print(final_df)

        print('------')
        print('adding nonsyn')

        print(final_df)
        print('--------')
        print('checking major minor')

        final_df['majRefSame'] = final_df.apply(lambda x: (x['major'] == x['refnt']), axis=1)
        final_df['minRefSame'] = final_df.apply(lambda x: (x['minor'] == x['refnt']), axis=1)
        #print(final_df)

        if "nonsyn" not in final_df: #the nonsyn column won't be added if no minor var
            final_df['nonsyn'] = ""
            final_df.to_csv("{0}/{1}.{2}.{3}.{4}.aa.snplist.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut),index=False)
            
        else:
            final_df.to_csv("{0}/{1}.{2}.{3}.{4}.aa.snplist.csv".format(args.save_dir,name,STRAIN,SEGMENT,freqcut),index=False)
