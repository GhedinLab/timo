## Written by: Kate Johnson, kej310@nyu.edu
"""
This code does three different things:
1. Checks the snplist file to make sure all positions are present, and not duplicated
        - if there is missing data (no reads aligning to position),
          then the position won't be present in snplist file
        - the updated csv files will be saved where the original snplist csv were

2. Generates a consensus sequence using checked/updated snplist file
-Generates by 'segment'
-If coverage is below provided cutoff, the nt will be replaced with an "N"

3. Pulls coverage information- that can then be used for coverage figures, etc

4. Will make separate fasta files for seqs

Input parameters:
--ref: path to reference, and name of reference fasta file

--var: path to the variant files (snplist)

--cov: coverage cutoff for calling a major nt. Default is 10x.

--minfreq: minorfreq cutoff used to call minor variants using the readreport script.
            This will be used to pull out files, and name files. So it must be
            the same as what was used to generate snplist files.

--strain: the strain named that was used to generate the snplist files. For example:
            COV19, H1N1, H3N2. These will also be present in the name of the snplist files.

-- savecov: The directory where you want to save your coverage files. Default is
            where the python script is running

--savecon: The directory where you want to save your consensus files. Default is
            where the python script is running.


TO RUN (examples):

python3 ConsensusFasta.Coverage.v3cov.py --ref ../MergedRuns/reference/SARS-COV-2.fasta
                                        --var ../MergedRuns/NYU/fullvarlist/
                                        --strain cov19
                                        --savecov ../MergedRuns/NYU/coverage
                                        --savecon ../MergedRuns/NYU/consensus


python3 ConsensusFasta.Coverage.v3cov.py --ref ../MergedRuns/reference/SARS-COV-2.fasta
                                        --var ../MergedRuns/NYU/fullvarlist/
                                        --strain cov19
                                        --savecov ../MergedRuns/NYU/coverage
                                        --savecon ../MergedRuns/NYU/consensus
                                        --cov 50
                                        --minfreq 0.02
"""


import os
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref','-r',required=True,help='Indicate path and reference fasta file')
parser.add_argument('--var','-v',required=True,help='Indicate path to variant files')
parser.add_argument('--cov','-c', default=10, help='indicate coverage cutoff threshold')
parser.add_argument('--minfreq','-f', default=0.001, help='indicate minfreq cutoff threshold')
parser.add_argument('--strain','-s',required=True,help='Indicate strain')
parser.add_argument('--savecov',default = '.', help='Indicate directory to save coverage csv')
parser.add_argument('--savecon',default = '.', help='Indicate directory to save coverage csv')
args = parser.parse_args()


def read_fasta(fp):
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
    segdict = {}
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            segdict[name[1:]] = seq
    return segdict


def add_row(df,ntpos_list,segment,name):
    """
    INPUT: Dataframe, and positions that differ from ref
    OUTPUT: Updated dataframe, with positions that were originally missing
    """
    for ntpos in ntpos_list: #iterate through list
        #generate new row of data
        new_row = {'name':name,'segment':segment,'ntpos':ntpos,'major':"N",'majorfreq':0,'minor':"",'minorfreq':"",'binocheck':"",'A':0,'C':0,'G':0,'T':0,'-':0,'totalcount':0,'aapos':"",'majoraa':"",'majorcodon':"",'minoraa':"",'minorcodon':""}
        #append to dataframe
        df = df.append(new_row,ignore_index=True)

    return df #return updated dataframe


#########################################################################################################################
ref = args.ref
path = args.var
STRAIN = args.strain
COVERAGE_CUTOFF = int(args.cov)#10 #the lowest coverage that you will allow - this is usually acceptable for the major nucleotide
refdict = open_fasta(ref)
freqcut = args.minfreq

cdf = pd.DataFrame()

for SEGMENT in refdict:
    print(SEGMENT)
    poor_alignments = []
    outf = open('{3}/{0}.{1}.{2}.fasta'.format(SEGMENT,COVERAGE_CUTOFF,STRAIN,args.savecon), 'w')
    outf.write('>' + SEGMENT + '\n' + refdict[SEGMENT] + '\n')
    seg_len = len(refdict[SEGMENT]) #length of segments
    #print(seg_len)

    #generate list so we pull positions that are within this range
    seg_list = list(range(1, seg_len + 1))

    #print(seg_list) #covid should be 29903
    #print(seg_list)
    for infile in glob.glob( os.path.join(path, '*{0}.{1}.{2}.snplist.csv'.format(STRAIN.upper(),SEGMENT, freqcut)) ):#will go through each specified csv file
        print(infile)
        filename = str(infile)
        filename_split = filename.split(SEGMENT)
        new_filename = "{0}{1}.Updated{2}".format(filename_split[0],SEGMENT,filename_split[1])


        df = pd.read_csv(infile) #pandas will read the csv file create dataframe

        df = df.drop_duplicates(subset=['ntpos'],keep='first') #drop any dups that are generated

        name = list(set(df['name']))[0] #pull name of sample
        print(name)


        #fullname = '{0}_{1}'.format(name, SEGMENT) #name of sequence for fasta file
        fullname = name

        df = df[df.ntpos.isin(seg_list)] #subset dataframe for ntpos only in region of interest

        ntpos_list = list(df['ntpos']) #pull nt pos from snplist file

        differ = list(set(seg_list) - set(ntpos_list)) #determine if there are regions missing
        print(sorted(differ))

        if len(differ) != 0 :
            #if the snplist file and reference differ in position, add those positions
            #add row using function add_row
            #df,ntpos_list,segment,name
            ndf = add_row(df,differ, SEGMENT, name)
            poor_alignments.append(name) #provide all the samples that need a new file

        else:
            #if there is no difference in snplist and reference then don't change
            ndf = df

        ndf = ndf.sort_values(by=['ntpos']) #sort df by ntpos

        ndf = ndf.reset_index(drop=True) #reset index and drop original index

        # write to an 'updated' snplist csv files
        ndf.to_csv(new_filename,index=False)

        # now moving on to pull coverage information
        # pull info that is useful for plotting coverage
        covdf = ndf[['name','segment','ntpos','totalcount']]

        # append information to empty df (usefull if there are actual segments)
        cdf = cdf.append(covdf)

        # generate consensus
        # change any low coverage positions to 'N' for the major
        ndf.loc[ndf.totalcount < COVERAGE_CUTOFF, "major"] = 'N'

        #generate consensus seq by pulling column in sorted,reindexed df
        consensus_seq = "".join(list(ndf['major']))

        print(len(consensus_seq))

        # output individual files
        outi = open('{3}/{0}.{1}.{2}.{4}.fasta'.format(name,COVERAGE_CUTOFF,STRAIN,args.savecon,SEGMENT), 'w')
        outi.write('>' + name + '\n' + consensus_seq + '\n')
        outi.close()

        #write to fasta file
        outf.write('>' + name + '\n' + consensus_seq + '\n')

    print(SEGMENT)

    print(poor_alignments)

    outf.close() #close fasta file (by segment)

# write the coverage data to dataframe
cdf.to_csv("{0}/{1}.coverage.csv".format(args.savecov,STRAIN.upper()),index=False)
