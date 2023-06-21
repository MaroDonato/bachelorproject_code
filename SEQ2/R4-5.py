'''
*** Running instructions ***

The script will look at a series of (one or more) FastQ files in a given folder, using the last character in the file name to organise them.
This last character should be a series of consecutive numbers. There should be no other FastQ files in that same folder.

1. Fill in below:
-the file path for the input file(s) as 'DataFolderLocation'
-the selection round that will be used for sequence sorting (number must also be at end of file name) as 'BaseSelectionRoundNumber'
-the number of sequences to extract per round as 'TopNPeptidesNumber'
-the output file name (and path) as 'SummaryFileName'
-your initiating amino acid

2. Run the code

3. Output files will be in the same folder as the input files

Notes:
-other reprogramming can be defined in the codon table further down under the function OpenReadingFrame; ATG elongator separate from initiator. All case sensitive one letter codes
-Assumes sequences have a GSGSGS* linker at the end (or similar length). If not, change definition of 'peptideParameters = TopNPeptidesList' before outputs and the primer sequence in 'SortedPeptideSequencesList'
-Proteases cut in the middle of the four letters given

Common problems:
-Make sure round number is present in file names
-Make sure input files are fastq format
-Make sure TopNPeptidesNumber is not higher than the number of viable sequences
-check file paths carefully, and remove any extra characters added when dragging and dropping (e.g. '\')
'''
import os
path = os.getcwd()
DataFolderLocation = r'/Users/maro/PycharmProjects/bachelorproject_sequences/SEQ2'
BaseSelectionRoundNumber = 4
TopNPeptidesNumber = 100
SummaryFileName = r'/Users/maro/PycharmProjects/bachelorproject_sequences/SEQ2/summary'
Initiator = 'Y'
expanded = 0    #change to 1 for expanded details about sequences, or 0 for simple sorting

#-------------------------------------------------------------------------------
# Creates a nested dictionary for each property in which each amino acid is linked to a value.
AminoAcidsDict = {
    'Name': 
    { 
      'A':'Alanine',        'R':'Arginine',         'N':'Aspargine',          'D':'Aspartic Acid',
      'C':'Cysteine',       'E':'Glutamic Acid',    'Q':'Glutamine',          'G':'Glycine',
      'H':'Histidine',      'O':'Hydroxy Proline',  'I':'Iso Leucine',        'L':'Leucine',
      'K':'Lysine',         'M':'Methionine',       'F':'Phenylalanine',      'P':'Proline',
      'U':'Pyroglutamic',   'S':'Serine',           'T':'Threonine',          'W':'Trypthophan',
      'Y':'Tyrosine',       'V':'Valine',           '*':'None'
    },

    'MolecularWeight':
    {
      'A':71.08,        'R':156.19,         'N':114.11,         'D':115.09,
      'C':103.15,       'E':219.12,         'Q':128.13,         'G':57.05,
      'H':137.14,       'O':113.11,         'I':113.16,         'L':113.16,
      'K':128.18,       'M':131.20,         'F':14.18,          'P':97.12,
      'U':121.09,       'S':87.08,          'T':101.11,         'W':186.22,
      'Y':163.22,       'V':99.13,          '*': 0.0
    },

    'pKx_SideChain': 
    {
      'A':None,         'R':12.48,          'N':None,           'D':3.56,
      'C':8.18,         'E':4.25,           'Q':None,           'G':None,
      'H':6.00,         'O':None,           'I':None,           'L':None,
      'K':10.53,        'M':None,           'F':None,           'P':None,
      'U':None,         'S':None,           'T':None,           'W':None,
      'Y':10.07,        'V':None,           '*':None
    },

    'ChargeP+':
    {
      'R':1,            'D':0,              'C':0,              'E':0,
      'H':1,            'K':1,              'Y':0
    },

    'ChargeP-':
    {
        'R':0,          'D':-1,             'C':-1,             'E':-1,
        'H':0,          'K':0,              'Y':-1
    },

    'GRAVYindex': 
    {
        '*':None,       'A':1.8,            'R':-4.5,           'N':-3.5,
        'D':-3.5,       'C':2.5,            'E':-3.5,           'Q':-3.5,
        'G':-0.4,       'H':-3.2,           'O':None,           'I':4.5,
        'L':3.8,        'K':-3.9,           'M':1.9,            'F':2.8,
        'P':-1.6,       'U':None,           'S':-0.8,           'T':-0.7,
        'W':-0.9,       'Y':-1.3,           'V':4.2
    },

    'Extinction_Coefficient280':
    {
      'W':5500,         'Y':1490,       'C':125
    },

    'Extinction_Coefficient214':
    {
        '*':None,           'A':32,         'R':102,            'N':136,
        'D':58,             'C':225,        'E':78,             'Q':142,
        'G':21,             'H':5125,       'O':0,              'I':45,
        'L':45,             'K':41,         'M':980,            'F':5200,
        'P':30,             'U':0,          'S':34,             'T':41,
        'W':29050,          'Y':5375,       'V':43
    },

    'PP1':
    {
        'A':-0.96,          'R':0.80,           'N':0.82,           'D':1.00,
        'C':-0.55,          'E':0.94,           'Q':0.78,           'G':-0.88,
        'H':0.67,           'O':0,              'I':-0.94,          'L':-0.90,
        'K':0.60,           'M':-0.82,          'F':-0.85,          'P':-0.81,
        'U':0,              'S':0.41,           'T':0.40,           'W':1.00,
        'Y':0.31,           'V':-1.00
    },

    'PP2':
    {
        'A':-0.76,          'R':0.63,           'N':-0.57,          'D':-0.89,
        'C':-0.47,          'E':0.54,           'Q':-0.30,          'G':-1.00,
        'H':-0.11,          'O':0,              'I':-0.05,          'L':0.03,
        'K':0.10,           'M':0.03,           'F':0.48,           'P':-0.40,
        'U':0,              'S':-0.82,          'T':-0.64,          'W':1.00,
        'Y':0.42,           'V':-0.43
    }
}

def CheckNetCharge(pH, PeptideLibraryList):
    """ Calculates the sum of the formal charge of each sequence at the given pH, by checking the pK of the side chain. Depending on its nature
    the amino acid is either protonated or deprotonated.
    """
    length = len(PeptideLibraryList)
    NetChargeSub = []
    NetChargeSum = []
    i = 0
    while i < length:
        for character in PeptideLibraryList[i]:
            if AminoAcidsDict['pKx_SideChain'][character] is None:
                continue
            elif pH < AminoAcidsDict['pKx_SideChain'][character]:  
                NetChargeSub.append(AminoAcidsDict['ChargeP+'][character])       
            elif pH > AminoAcidsDict['pKx_SideChain'][character]:
                NetChargeSub.append(AminoAcidsDict['ChargeP-'][character])
            elif pH == AminoAcidsDict['pKx_SideChain'][character]:
                NetChargeSub.append((AminoAcidsDict['ChargeP+'][character] + AminoAcidsDict['ChargeP-'][character]) / 2)
        NetChargeSub.append(None) ## none used as seperator of each sequence. coulnd figure out how to append each sequences charge to a sublist right away
        i += 1

    arrays = [[]]
    for i, item in enumerate(NetChargeSub):
        arrays[-1].append(item)
        if item == None and i != len(NetChargeSub)-1:
            arrays.append([])

    NetChargeClean = []
    NetChargeClean = [[i for i in nested if i != None] for nested in arrays]
    ## print(len(NetChargeClean) == len(arrays)) 


    max_len = max(len(item) for item in NetChargeClean) # to sum the charges to get the net charge each nested list in NetchargeClean has to have the same amount of items.

    for item in NetChargeClean:
        if len(item) < max_len:
            item.extend([0] * (max_len - len(item))) # Extends each nested list with 0 to equal the longest list.

    NetChargeSum = [sum(x) for x in NetChargeClean]
    #print (NetChargeSum)
    return NetChargeSum

""" COUNT THE LOOPS AND ITS SIZE """
def Loops_Size(PeptideLibraryList):
    """ Return a nested list ind which [A,B,C] A= the amount of loops in the amino acid 
    B= is the size of the first loop and C= is the size of the second loop"""

    C_amount = []
    Size_Loop1 = []
    Size_Loop2 = []

    for _ in range(len(PeptideLibraryList)):
        Size_Loop1.append(None)
    for _ in range(len(PeptideLibraryList)):
        Size_Loop2.append(None)


    for item in PeptideLibraryList:
        C_amount.append((item.count("C", 2, len(item)))) #starts counting at the index of 2 which skips the first and second amino acid as it is not able to create loops

    for ind, value in enumerate(C_amount):
        if value == 1:
            s = PeptideLibraryList[ind]
            c = s.index("C") 
            Size_Loop1.pop(ind)
            Size_Loop1.insert(ind, c)
        elif value == 2:
            s = PeptideLibraryList[ind]
            c = s.index("C") 
            Size_Loop1.pop(ind)
            Size_Loop1.insert(ind, c)
        elif value == 3:
            s = PeptideLibraryList[ind]
            c = s.index("C") #first C
            c2 = s.index("C", c+1) #second C
            c3 = s.index("C", c2+1)
            size1 = c
            size2 = c3 - c2
            Size_Loop1.pop(ind)
            Size_Loop1.insert(ind, size1)
            Size_Loop2.pop(ind)
            Size_Loop2.insert(ind, size2)
        elif value == 4:
            s = PeptideLibraryList[ind]
            c = s.index("C") #first C
            c2 = s.index("C", c+1) #second C
            c3 = s.index("C", c2+1)
            size1 = c
            size2 = c3 - c2
            Size_Loop1.pop(ind)
            Size_Loop1.insert(ind, size1)
            Size_Loop2.pop(ind)
            Size_Loop2.insert(ind, size2)
        elif value > 4:
            Size_Loop1.pop(ind)
            Size_Loop1.insert(ind, 'more than 2 loops observed') #how far? is relevant?

    Amount_Loops = []
    for _ in range(len(PeptideLibraryList)):
        Amount_Loops.append(None)

    for ind, value in enumerate(C_amount):
        if value == 1:
            Amount_Loops.pop(ind)
            Amount_Loops.insert(ind, 1)
        else:
            a = (value + 1) // 2
            Amount_Loops.pop(ind)
            Amount_Loops.insert(ind, a)
    
    Loops_Size= []
    Loops_Size= [list(e) for e in zip(Amount_Loops, Size_Loop1, Size_Loop2)]    
    #FinalLoops = [" ".join([str(c) for c in lst]) for lst in Loops_Size]
    
    return Loops_Size

def AmountLip(PeptideLibraryList):
    """"Counts the total amount of amino acids and returns the amount of which are lipophilic as a ration x out of total"""
    import re
    Amount_Lipophilic = []
    List_matches = []
    for item in PeptideLibraryList:
        l = len(item)
        pattern = re.compile('[AVLIPFMWCY]') #Hydrophobic aminoacids
        List_matches = pattern.findall(item)
        count = len(List_matches)
        Amount_Lipophilic.append(str(count) + " out of "+str(l))
    return Amount_Lipophilic

""" GRAVY HYDROPATHY """

def GRAVY(PeptideLibraryList):
    """Caculated the GRAVY index ref:
    [1] Kyte J, Doolittle RF, Diego S, Jolla L. A Simple Method for Displaying the Hydropathic Character of a Protein 1982:105–32. """
    GRAVY = []
    for item in PeptideLibraryList:
        Values = []
        l = len(item)
        for character in item:
            if character == '*':
                pass
            else:    
                Values.append(AminoAcidsDict['GRAVYindex'][character])
        r = sum(Values)/l
        GRAVY.append(r)
        Values.clear()
    
    GRAVY = ['%.2f' % elem for elem in GRAVY]
    return GRAVY

""" MOLAR EXTINCTION COEFFICIENT AT 280 nm """
def E_280(PeptideLibraryList):
    import re
    E_280 = []
    Matches = []
    for _ in range(len(PeptideLibraryList)):
        E_280.append(0)
    for index, item in enumerate(PeptideLibraryList):
        pattern = re.compile('[WYC]') #Hydrophobic aminoacids
        Matches = pattern.findall(item)
        W = Matches.count('W')
        Y = Matches.count('Y')
        C = (Matches.count('C') - 1)//2 # only adds a value of 250 for a c-c pair. One cysteine is skipped as it already is used in the cyclysation.
        
        if C >= 0:
            Sum = (W*5500)+(Y*1490)+(C*250)
        else:
            Sum = (W*5500)+(Y*1490)
        
        del E_280[index]
        E_280.insert(index, Sum)

    return E_280

def E_214(PeptideLibraryList):
    E_214 = []
    for item in PeptideLibraryList:
        Values = []
        for character in item:
            if character == '*':
                pass
            else:
                Values.append(AminoAcidsDict['Extinction_Coefficient214'][character])
        E_214.append(sum(Values))
        Values.clear()
    return E_214

def CellPermiability(PeptideLibraryList):
    """ Return a list with item = each peptide and either falls within the green area = yes, yellow area = maybe or red are = no
    If the function returns Yes and or maybe it also gives the values voor (x.y) = (PP2, PP1) = (Hydrophobicity, Polarity)
    Reference: [1] Schmidt S, Adjobo-hermans MJW, Kohze R, Enderle T, Brock R, Milletti F. Identi fi cation of Short Hydrophobic Cell-Penetrating Peptides for Cytosolic Peptide Delivery by Rational Design 2017. doi:10.1021/acs.bioconjchem.6b00535."""


    PP1 = []
    PP2 = []
    temp2 = []
    temp1 = []
    for item in PeptideLibraryList:
        for character in item:
            if character == '*':
                pass
            else:
                temp1.append(AminoAcidsDict["PP1"][character])
                temp2.append(AminoAcidsDict['PP2'][character])
        
        PP1.append(sum(temp1)/len(item))
        PP2.append(sum(temp2)/len(item))

        temp1=[]
        temp2=[]
    PP1= ['%.2f' % elem for elem in    PP1]
    PP1 = [float(elem) for elem in PP1]
    PP2= ['%.2f' % elem for elem in PP2]
    PP2 = [float(elem) for elem in PP2]
    
    # (x,y) 
    PP2_PP1 = zip(PP2,   PP1)
    PP2_PP1 = list(PP2_PP1)
    

    CellPermiability = []
    for item in PP2_PP1:
        if item[0] >= 0:
            CellPermiability.append("Yes, "+ str(item))
        if item[0] < 0:
            if item[1] <= (1.5*item[0]-0.25):
                CellPermiability.append("Yes, " + str(item))
            elif item[1] <= (1.5*item[0]+0.25):
                CellPermiability.append("Maybe, " + str(item))
            else:
                CellPermiability.append("No,,")

    #FinalPerm = [" ".join([str(c) for c in lst]) for lst in CellPermiability]

    return CellPermiability
            

"""GIVE MOLECULAR WEIGHT """
def MolecularWeight(peptideSequence):

    #PeptideLibraryListMWall = []
    #i = 0
    #while i < length:
    peptidemass = 41.05+16.02 #assuming N-terminal chloroacetyl group cyclised onto Cys and C-terminal acetamide
    numCys = peptideSequence.count("C", 2, len(peptideSequence))
    if numCys == 3:
        peptidemass -= 2 #accounting for disulfide formation
    for character in peptideSequence:
        if character == "*":
            break
        else:
            peptidemass+=AminoAcidsDict['MolecularWeight'][character] 
     #   i += 1 ## Translates the Peptide code to its MW and puts it in a List

    #PeptideMW = [ sum(PeptideLibraryListMWall[x:x+23]) for x in range(0, len(PeptideLibraryListMWall), 23)] ## 23 is the length of each sequence, more elegant way to devide eacah sequence?
    #PeptideMWrounded = [ '%.2f' % elem for elem in peptidemass] # gives two point decimals of the previous list
    return(peptidemass)

def peptideLength(peptideSequence):
    '''
    Takes a peptide sequence and returns its length, up to the final Cys
    '''
    outputLength = len(peptideSequence)-7 #assumes GSGSGS* at end of each sequence, could be tested for and handled more elegantly for truncations in linker
    return (outputLength)

def pI(PeptideLibraryList):
    """returns ranges of pH in which the peptides charge is 0, a bit faster than the previous function"""
    import numpy
    pHrange = []
    for x in numpy.arange(0.0, 14.1, 0.01):
        pHrange.append(x)
    pHrange = ['%.2f' % elem for elem in pHrange]
    pHrange = [float(elem) for elem in pHrange]
    
    def CheckNetChargeString(string, pH):
        
        NetChargeSub = []
        for character in string:
            
            try:
                if AminoAcidsDict['pKx_SideChain'][character] == None:
                    pass
                elif pH < AminoAcidsDict['pKx_SideChain'][character]:  
                    NetChargeSub.append(AminoAcidsDict['ChargeP+'][character])       
                elif pH > AminoAcidsDict['pKx_SideChain'][character]:
                    NetChargeSub.append(AminoAcidsDict['ChargeP-'][character])
                elif pH == AminoAcidsDict['pKx_SideChain'][character]:
                    NetChargeSub.append((AminoAcidsDict['ChargeP+'][character] + AminoAcidsDict['ChargeP-'][character]) / 2)
            except KeyError:
                pass

        NetChargeSum = sum(NetChargeSub)
        
        return NetChargeSum

    pI = []
    temp = []
    for item in PeptideLibraryList:
    
        for value in pHrange:
            if CheckNetChargeString(item, value) == 0:
                temp.append(value)
            else:
                pass
        if len(temp) > 1:
            pI.append([temp[0], temp[-1]])
        else:
            pI.append(['None', 'None'])
            
        temp = []

    #FinalpI = [" ".join([str(c) for c in lst]) for lst in pI]

    return pI


def Trypsin(PeptideLibraryList):      
    """ Returns a list per peptide on the sites where the peptidase might act upon"""

    Trypsin = [] 
    per_itemlist = []
    Exceptions_trypsin=['CKD', 'DKD', 'CKH', 'CKY','CRK','RRH','RRR']  
    for item in PeptideLibraryList:
        try:
            K_count = item.count('K',1)
            R_count = item.count('R',1)

            if K_count + R_count == 0:
                Trypsin.append("None")
            else:
                a = 1
                a = item.find('K', a, -1)
                b = 1
                b = item.find('R', b, -1)
                while K_count > 0:
                    e1 = item[a-1:a+2] 
                    if e1 in Exceptions_trypsin:
                        pass              
                    elif (item[a+1] == 'P') and (item[a-1] == 'W'):
                            per_itemlist.append(item[a-1:a+3])   
                    else:
                        if item[a+2] != "P":
                            per_itemlist.append(item[a-1:a+3])
            
                    a = item.find('K', a+1 , -1)
                    K_count = K_count - 1

                while R_count > 0:
                    e2 = item[b-1:b+2] 
                    if e2 in Exceptions_trypsin:
                        pass              
                    elif (item[b+1] == 'P') and (item[b-1] == 'M'):
                            per_itemlist.append(item[b-1:b+3])   
                    else:
                        if item[b+2] != "P":
                            per_itemlist.append(item[b-1:b+3])
            
                    b = item.find('R', b+1 , -1)
                    R_count = R_count - 1
                if len(per_itemlist) > 0:
                    Trypsin.append(per_itemlist)
                    per_itemlist = []
                else:
                    Trypsin.append("None")
        except IndexError:
            Trypsin.append("error")

    FinalTrypsin = [" ".join([str(c) for c in lst]) for lst in Trypsin]

    return FinalTrypsin

def Pepsin_pH1_3(PeptideLibraryList):
    """ Pepsin at pH 1.3: Returns a list per peptide on the sites where the peptidase might act upon""" 

    Pepsin = []
    per_itemlist = []
    exceptions_Pepsin = ['HPR', 'KPR', 'RPR', 'HP', 'KP', 'RP']
    for item in PeptideLibraryList:
        try:
            F_count = item.count('F',1)
            L_count = item.count('L',1)

            if F_count + L_count == 0:
                Pepsin.append("None")
            else:
                a = 1
                a = item.find('F', a, -1)
                b = 1
                b = item.find('L', b, -1)
                while F_count > 0:
                    if (item[a+1] == 'P') or (item[a+2]  == 'P'):
                        pass  
                    elif item[a-2:a] not in exceptions_Pepsin:
                        per_itemlist.append(item[a-1:a+3])
                    elif item[a-3:a] not in exceptions_Pepsin:
                        per_itemlist.append(item[a-2:a+2])
                    else:
                        continue
                    a = item.find('F', a+1)
                    F_count = F_count - 1

                while L_count > 0:
                    if (item[b+1] == 'P') or (item[b+2]  == 'P'):
                        pass  
                    elif item[b-2:b] not in exceptions_Pepsin:
                        per_itemlist.append(item[b-1:b+3])
                    elif item[b-3:b] not in exceptions_Pepsin:
                        per_itemlist.append(item[b-2:b+2])
                    else:
                        continue
            
                    b = item.find('L', b+1 , -1)
                    L_count = L_count - 1
                if len(per_itemlist) > 0:
                    Pepsin.append(per_itemlist)
                    per_itemlist = []
                else:
                    Pepsin.append("None")
        except IndexError:
            Pepsin.append("error")
    FinalPepsin = [" ".join([str(c) for c in lst]) for lst in Pepsin]

    return FinalPepsin



def Pepsin_pH3(PeptideLibraryList):
    """Pepsin pH > 3: Returns a list per peptide on the sites where the peptidase might act upon """

    Pepsin = []
    per_itemlist = []
    exceptions_Pepsin = ['HPR', 'KPR', 'RPR', 'HP', 'KP', 'RP']
    for item in PeptideLibraryList:
        try:
            F_count = item.count('F',1)
            L_count = item.count('L',1)
            W_count = item.count('W',1)
            Y_count = item.count('Y',1)

            if F_count + L_count + W_count + Y_count == 0:
                Pepsin.append("None")
            else:
                a = 1
                a = item.find('F', a, -1)
                b = 1
                b = item.find('L', b, -1)
                c = 1
                c = item.find('W', c, -1)
                d = 1
                d = item.find('Y', d, -1)
                while F_count > 0:
                    if (item[a+1] == 'P') or (item[a+2]  == 'P'):
                        pass  
                    elif item[a-2:a] not in exceptions_Pepsin:
                        per_itemlist.append(item[a-1:a+3])
                    elif item[a-3:a] not in exceptions_Pepsin:
                        per_itemlist.append(item[a-2:a+2])
                    else:
                        continue
                    a = item.find('F', a+1)
                    F_count = F_count - 1

                while L_count > 0:
                    if (item[b+1] == 'P') or (item[b+2]  == 'P'):
                        pass  
                    elif item[b-2:b] not in exceptions_Pepsin:
                        per_itemlist.append(item[b-1:b+3])
                    elif item[b-3:b] not in exceptions_Pepsin:
                        per_itemlist.append(item[b-2:b+2])
                    else:
                        continue
            
                    b = item.find('L', b+1 , -1)
                    L_count = L_count - 1

                while W_count > 0:
                    if (item[c+1] == 'P') or (item[c+2]  == 'P'):
                        pass  
                    elif item[c-2:c] not in exceptions_Pepsin:
                        per_itemlist.append(item[c-1:c+3])
                    elif item[c-3:c] not in exceptions_Pepsin:
                        per_itemlist.append(item[c-2:c+2])
                    else:
                        continue
            
                    c = item.find('W', c+1)
                    W_count = W_count - 1
                
                while Y_count > 0:
                    if (item[d+1] == 'P') or (item[d+2]  == 'P'):
                        pass  
                    elif item[d-2:d] not in exceptions_Pepsin:
                        per_itemlist.append(item[d-1:d+3])
                    elif item[d-3:d] not in exceptions_Pepsin:
                        per_itemlist.append(item[d-2:d+2])
                    else:
                        continue
            
                    d = item.find('Y', d+1)
                    Y_count = Y_count - 1

                if len(per_itemlist) > 0:
                    Pepsin.append(per_itemlist)
                    per_itemlist = []
                else:
                    Pepsin.append("None")
        except IndexError:
            Pepsin.append("error")
        
    FinalPepsin = [" ".join([str(c) for c in lst]) for lst in Pepsin]

    return FinalPepsin
    
def Chymotrypsin_HighSpecificity(PeptideLibraryList):
    Chemotrypsin_HighSpecificity = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            p = ['F','Y','W']
            matches = [x for x in p if x in item]
            p_count =len(matches)
            if p_count == 0:
                Chemotrypsin_HighSpecificity.append("None")
            else:
                for value in p:
                    if value == 'W':
                        a = 1
                        a = item.find(value, a, -1)
                        p_count = item.count(value, 1)
                        while p_count > 0:
                            if a > 1:
                                if (item[a+1] != 'M') and (item[a+1] != 'P'):
                                    per_itemlist.append(item[a-1:a+3])
                                else:
                                    pass
                            else:
                                pass
                            a = item.find(value, a+1,-1)
                            p_count -= 1
                    else:
                        a = 1
                        a = item.find(value, a, -1)
                        p_count = item.count(value, 1)
                        while p_count > 0:
                            if a > 1:
                                if (item[a+1] != 'P'):
                                    per_itemlist.append(item[a-1:a+3])
                                else:
                                    pass
                            else:
                                pass
                            a = item.find(value, a+1,-1)
                            p_count -= 1

                if len(per_itemlist) > 0:
                    Chemotrypsin_HighSpecificity.append(per_itemlist)
                    per_itemlist = []
                else:
                    Chemotrypsin_HighSpecificity.append("None")
        except IndexError:
            Chemotrypsin_HighSpecificity.append("error")
    FinalChymo = [" ".join([str(c) for c in lst]) for lst in Chemotrypsin_HighSpecificity]

    return FinalChymo

def Chymotrypsin_LowSpecificity(PeptideLibraryList):
    """ M is not considered in this function"""

    Chemotrypsin_LowSpecificity = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            p = ['F','L','Y','W','H']
            patternH = ['D','M','P','W']
            matches = [x for x in p if x in item]
            p_count =len(matches)
            if p_count == 0:
                Chemotrypsin_LowSpecificity.append("Peptide is not hydrolysed by Chemitrypsin Lowspecifity")
            else:
                for value in p:
                    if value == 'H':
                        a = 1
                        a = item.find(value, a, -1)
                        h_count = item.count(value, 1)
                        while h_count > 0:
                            if a > 1:
                                if item[a+1] not in patternH:
                                    per_itemlist.append(item[a-1:a+3])
                                else:
                                    pass
                            else:
                                pass
                            a = item.find(value, a+1,-1)
                            h_count -= 1
                    elif value == 'W':
                        a = 1
                        a = item.find(value, a, -1)
                        w_count = item.count(value, 1)
                        while w_count > 0:
                            if a > 1:
                                if (item[a+1] != 'P') or (item[a+1] != 'M'):
                                    per_itemlist.append(item[a-1:a+3])
                                else:
                                    pass
                            else:
                                pass
                            a = item.find(value, a+1,-1)
                            w_count -= 1
                    elif value == 'M':
                        a = 1
                        a = item.find(value, a, -1)
                        m_count = item.count(value, 1)
                        while m_count > 0:
                            if a > 1:
                                if (item[a+1] != 'P') or (item[a+1] != 'Y'):
                                    per_itemlist.append(item[a-1:a+3])
                                else:
                                    pass
                            else:
                                pass
                            a = item.find(value, a+1,-1)
                            m_count -= 1
                    else:
                        a = 1
                        a = item.find(value, a, -1)
                        p_count = item.count(value, 1)
                        while p_count > 0:
                            if a > 1:
                                if (item[a+1] != 'P'):
                                    per_itemlist.append(item[a-1:a+3])
                                else:
                                    pass
                            else:
                                pass
                            a = item.find(value, a+1,-1)
                            p_count -= 1

                if len(per_itemlist) > 0:
                    Chemotrypsin_LowSpecificity.append(per_itemlist)
                    per_itemlist = []
                else:
                    Chemotrypsin_LowSpecificity.append("None")
        except IndexError:
            Chemotrypsin_LowSpecificity.append("error")

    FinalChymo = [" ".join([str(c) for c in lst]) for lst in Chemotrypsin_LowSpecificity]

    return FinalChymo

def Misc_peptidases(PeptideLibraryList):
    """ Creates a list of every other peptidase on every item. """
    import re
    

    #Arg-C proteinase
    Arg_C = []
    per_itemlist =[]
    for item in PeptideLibraryList:
        try:
            R_count = item.count('R',1)
            if R_count == 0:
                Arg_C.append(['XXXX'])
            else:
                a = 1
                a = item.find('R', a, -1)
                while R_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('R',a+1,-1)
                    R_count -= 1
                if len(per_itemlist) > 0:
                    Arg_C.append(per_itemlist)
                    per_itemlist = []
                else:
                    Arg_C.append(['XXXX'])
        except IndexError:
            Arg_C.append(['XXXX'])
            
            
    
    # Asp-N endopeptidase
    Asp_N = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            D_count = item.count('D',1)
            if D_count == 0:
                Asp_N.append(['XXXX'])
            else:
                a = 1
                a = item.find('D', a, -1)
                while D_count > 0:
                    per_itemlist.append(item[a-2:a+2])
                    a = item.find('D',a+1,-1)
                    D_count -= 1

                if len(per_itemlist) > 0:
                    Asp_N.append(per_itemlist)
                    per_itemlist = []
                else:
                    Asp_N.append(['XXXX'])
        except IndexError:
            Asp_N.append(['XXXX'])

    # BNPS-Skatole
    BNPS = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            W_count = item.count('W',1)
            if W_count == 0:
                BNPS.append(['XXXX'])
            else:
                a = 1
                a = item.find('W', a, -1)
                while W_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('W',a+1,-1)
                    W_count -= 1
                if len(per_itemlist) > 0:
                    BNPS.append(per_itemlist)
                    per_itemlist = []
                else:
                    BNPS.append(['XXXX'])
        except IndexError:
            BNPS.append(['XXXX'])
    
    # Clostridiopeptidase B
    Clost = []
    per_itemlist =[]
    for item in PeptideLibraryList:
        try:
            R_count = item.count('R',1)
            if R_count == 0:
                Clost.append(['XXXX'])
            else:
                a = 1
                a = item.find('R', a, -1)
                while R_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('R',a+1,-1)
                    R_count -= 1
                if len(per_itemlist) > 0:
                    Clost.append(per_itemlist)
                    per_itemlist = []
                else:
                    Clost.append(['XXXX'])
        except IndexError:
            Clost.append(['XXXX'])

    # CNBr
    CNBr= []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            M_count = item.count('M',1)
            if M_count == 0:
                CNBr.append(['XXXX'])
            else:
                a = 1
                a = item.find('M', a, -1)
                while M_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('M',a+1,-1)
                    M_count -= 1
                if len(per_itemlist) > 0:
                    CNBr.append(per_itemlist)
                    per_itemlist = []
                else:
                    CNBr.append(['XXXX'])
        except IndexError:
            CNBr.append(['XXXX'])
    
    # Formic Acid
    Formic_Acid = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            D_count = item.count('D',1)
            if D_count == 0:
                Formic_Acid.append(['XXXX'])
            else:
                a = 1
                a = item.find('D', a, -1)
                while D_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('D',a+1,-1)
                    D_count -= 1
                if len(per_itemlist) > 0:
                    Formic_Acid.append(per_itemlist)
                    per_itemlist = []
                else:
                    Formic_Acid.append(['XXXX'])
        except IndexError:
            Formic_Acid.append(['XXXX'])
    
    # Glutamylendopeptidase
    Glutamyl = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            E_count = item.count('E',1)
            
            if E_count == 0:
                Glutamyl.append(['XXXX'])
            else:
                a = 1
                a = item.find('E', a, -1)
                while E_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('E',a+1,-1)
                    E_count -= 1
                if len(per_itemlist) > 0:
                    Glutamyl.append(per_itemlist)
                    per_itemlist = []
                else:
                    Glutamyl.append(['XXXX'])
        except IndexError:
            Glutamyl.append(['XXXX'])
    
    # Iodoso benzoic acid
    Iodoso = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            W_count = item.count('W',1)
            if W_count == 0:
                Iodoso.append(['XXXX'])
            else:
                a = 1
                a = item.find('W', a, -1)
                while W_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('W',a+1,-1)
                    W_count -= 1
                if len(per_itemlist) > 0:
                    Iodoso.append(per_itemlist)
                    per_itemlist = []
                else:
                    Iodoso.append(['XXXX']) 
        except IndexError:
            Iodoso.append(['XXXX']) 

    # LysC
    LysC = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            K_count = item.count('K',1)
            if K_count == 0:
                LysC.append(['XXXX'])
            else:
                a = 1
                a = item.find('K', a, -1)
                while K_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('K',a+1,-1)
                    K_count -= 1
                if len(per_itemlist) > 0:
                    LysC.append(per_itemlist)
                    per_itemlist = []
                else:
                    LysC.append(['XXXX'])
        except IndexError:
            LysC.append(['XXXX'])

    # NTCB (2-nitro-5thiocyanobezoic acid)
    NTCB = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            C_count = item.count('C',1)
            if C_count == 0:
                NTCB.append(['XXXX'])
            else:
                a = 1
                a = item.find('C', a, -1)
                while C_count > 0:
                    per_itemlist.append(item[a-2:a+2])
                    a = item.find('C',a+1,-1)
                    C_count -= 1
                if len(per_itemlist) > 0:
                    NTCB.append(per_itemlist)
                    per_itemlist = []
                else:
                    NTCB.append(['XXXX'])
        except IndexError:
            NTCB.append(['XXXX'])
    
    # Enterokinase
    
    Enterokinase = []
    pattern = ['DDD','EEE','DEE', 'EDE','EED', 'DDE','DED','EDD',]
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            K_count = item.count('K',1)
            if K_count == 0:
                Enterokinase.append(['XXXX'])
            else:
                per_itemlist = []
                a = 1
                a = item.find('K', a, -1)
                while K_count > 0:
                    if item[a-3:a] in pattern:
                        per_itemlist.append(item[a-1:a+3])
                    else:
                        pass
                    a = item.find('K',a+1,-1)
                    K_count -= 1
            
                if len(per_itemlist) > 0:
                    Enterokinase.append(per_itemlist)
                    per_itemlist = []
                else:
                    Enterokinase.append(['XXXX'])
        except IndexError:
            Enterokinase.append(['XXXX'])
        
    # Factor Xa
    FactorXa = []
    pattern = ['A','F','G','I','L','T','V','M']
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            R_count = item.count('GR',1)
            if R_count == 0:
                FactorXa.append(['XXXX'])
            else:
                a = 1
                a = item.find('GR', a, -1)
                while R_count > 0:
                    if (item[a-1] == 'D') or (item[a-1] == 'E'):
                        if item[a-2] in pattern:
                            per_itemlist.append(item[a:a+4])
                        else:
                            pass
                    else:
                        pass
                    a = item.find('GR',a+1,-1)
                    R_count = R_count - 1

                if len(per_itemlist) > 0:
                    FactorXa.append(per_itemlist)
                    per_itemlist = []
                else:
                    FactorXa.append(['XXXX'])
        except IndexError:
            FactorXa.append(['XXXX'])

    # Granzyme B
    GranzymeB = []
    per_itemlist= []
    for item in PeptideLibraryList:
        try:
            C_count = item.count('IEPD',1)
            if C_count == 0:
                GranzymeB.append(['XXXX'])
            else:
                a = 1
                a = item.find('IEPD', a, -1)
                while C_count > 0:
                    per_itemlist.append(item[a+2:a+6])
                    a = item.find('C',a+1,-1)
                    C_count -= 1

                if len(per_itemlist) > 0:
                    GranzymeB.append(per_itemlist)
                    per_itemlist = []
                else:
                    GranzymeB.append(['XXXX'])
        except IndexError:
            GranzymeB.append(['XXXX'])

    #Neutrophil elastase
    Elastase = []
    per_itemlist= []
    for item in PeptideLibraryList:
        try:
            A_count = item.count('A',1)
            V_count = item.count('V',1)
            if A_count + V_count == 0:
                Elastase.append(['XXXX'])
            else:
                a = 1
                a = item.find('A', a, -1)
                b = 1
                b = item.find('V', b, -1)
                while A_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('A',a+1,-1)
                    A_count -= 1
                
                while V_count > 0:
                    per_itemlist.append(item[b-1:b+3])
                    b = item.find('V',b+1,-1)
                    V_count -= 1

                if len(per_itemlist) > 0:
                        Elastase.append(per_itemlist)
                        per_itemlist = []
                else:
                    Elastase.append(['XXXX'])
        except IndexError:
            Elastase.append(['XXXX'])
    
    # Hydroxylamine
    Hydroxylamine = []
    per_itemlist= []
    for item in PeptideLibraryList:
        try:
            NG_count = item.count('NG',1)
            if NG_count == 0:
                Hydroxylamine.append(['XXXX'])
            else:
                a = 1
                a = item.find('NG', a, -1)
                while NG_count > 0:
                    per_itemlist.append(item[a-1:a+3])
                    a = item.find('NG',a+1,-1)
                    NG_count -= 1

                if len(per_itemlist) > 0:
                    Hydroxylamine.append(per_itemlist)
                    per_itemlist = []
                else:
                    Hydroxylamine.append(['XXXX'])
        except IndexError:
            Hydroxylamine.append(['XXXX'])
    
    
    #Proline Endopeptidase
    Proline_Endopeptidase = []
    import re
    pattern = re.compile('[HKR]')
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            P_count = item.count('P',1)
            if P_count == 0:
                Proline_Endopeptidase.append(['XXXX'])
            else:
                a = 1
                a = item.find('P', a, -1)
                while P_count > 0:
                    if (item[a+1] != 'P'):
                        b = item[a-1]
                        c = pattern.findall(b)
                        if len(c) > 0:
                            per_itemlist.append(item[a-1:a+3])
                        else:
                            pass
                    a = item.find('P',a+1,-1)
                    P_count = P_count - 1

                if len(per_itemlist) > 0:
                    Proline_Endopeptidase.append(per_itemlist)
                    per_itemlist = []
                else:
                    Proline_Endopeptidase.append(['XXXX'])
        except IndexError:
            Proline_Endopeptidase.append(['XXXX'])

    # Proteinase K
    ProteinaseK = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            p = ['A','E','F','I','L','V','W','Y']
            matches = [x for x in p if x in item]
            p_count = len(matches)
            if p_count == 0:
                ProteinaseK.append(['XXXX'])
            else:
                for value in p:
                    a = 1
                    a = item.find(value, a, -1)
                    p_count = item.count(value, 1)
                    while p_count > 0:
                        per_itemlist.append(item[a-1:a+3])
                        a = item.find(value, a+1,-1)
                        p_count -= 1

                if len(per_itemlist) > 0:
                    ProteinaseK.append(per_itemlist)
                    per_itemlist = []
                else:
                    ProteinaseK.append(['XXXX'])
        except IndexError:
            ProteinaseK.append(['XXXX'])
        
    

    #Thermolysin
    Thermolysin = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            p = ['A','F','I','L','V','M']
            matches = [x for x in p if x in item]
            p_count =len(matches)
            if p_count == 0:
                Thermolysin.append(['XXXX'])
            else:
                for value in p:
                    a = 1
                    a = item.find(value, a, -1)
                    p_count = item.count(value, 1)
                    while p_count > 0:
                        if a > 2:
                            if (item[a-1] != 'D') and (item[a-1] != 'E'):
                                per_itemlist.append(item[a-2:a+2])
                            else:
                                pass
                        else:
                            pass
                        a = item.find(value, a+1,-1)
                        p_count -= 1

                if len(per_itemlist) > 0:
                    Thermolysin.append(per_itemlist)
                    per_itemlist = []
                else:
                    Thermolysin.append(['XXXX'])
        except IndexError:
            Thermolysin.append(['XXXX'])
    
    # Staphylococcal peptidase I
    Staph = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            E_count = item.count('E',1)
            if E_count == 0:
                Staph.append(['XXXX'])
            else:
                a = 1
                while E_count > 0:              
                    a = item.find('E', a, -1)
                    if (item[a-1] != 'E'):
                        per_itemlist.append(item[a-1:a+3])
                    else:
                        pass
                    a = item.find('E',a+1,-1)
                    E_count = E_count - 1

                if len(per_itemlist) > 0:
                    Staph.append(per_itemlist)
                    per_itemlist = []
                else:
                    Staph.append(['XXXX'])
        except IndexError:
            Staph.append(['XXXX'])
    
    # Thrombin
    Thrombin = []
    per_itemlist = []
    pattern1 = ['A','F','G','I','L','T','V','W','A']
    pattern2 = ['A','F','G','I','L','T','V', 'M']
    for item in PeptideLibraryList:
        try:
            PR_count = item.count('PR',2)
            GRG_count = item.count('GRG')
            if PR_count + GRG_count == 0:
                Thrombin.append(['XXXX'])
            else:
                a = 1
                while GRG_count > 0:   
                    a = item.find('GRG',a ,-1)
                    per_itemlist.append(item[a:a+4])
                    a = item.find('GRG', a+1, -1)
                    GRG_count = GRG_count - 1
                a = 2
                while PR_count > 0:
                    
                    a = item.find('PR', a, -1)
                    if item[a-1] in pattern1:
                        if item[a-1] in pattern2:
                            if (item[a+2] == 'E') or (item[a+2] == 'D'):
                                pass
                            elif (item[a+3] == 'E') or (item[a+3] == 'D'):
                                pass
                            else:
                                per_itemlist.append(item[a:a+4])
                        else:
                            pass
                    else:
                        pass
                    a = item.find('PR',a+1,-1)
                    PR_count = PR_count - 1

                if len(per_itemlist) > 0:
                    Thrombin.append(per_itemlist)
                    per_itemlist = []
                else:
                    Thrombin.append(['XXXX'])
        except IndexError:
            Thrombin.append(['XXXX'])
    
    #Caspase1

    caspase1 = []
    per_itemlist =[]
    pattern =['H','A', 'T']
    pattern2 = ['P','E','D','Q','K','R']
    pattern3 = ['F','W','Y','L']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('D',3)
            if x_count == 0:
                caspase1.append(['XXXX'])
            else:            
                a = 3
                while x_count > 0:
                    a = item.find('D', a, -1)
                    if item[a-1] in pattern:
                        if item[a-3] in pattern3:
                            if item[a+1] not in pattern2:
                                per_itemlist.append(item[a-1:a+3])
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                    a = item.find('D',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase1.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase1.append(['XXXX'])
        except IndexError:
            caspase1.append(['XXXX'])
    
  

    # Caspase 2
    caspase2 = []
    per_itemlist = []
    pattern = ['P','E','D','Q','K','R']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('DVAD',1)
            if x_count == 0:
                caspase2.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('DVAD', a, -1)           
                    if item[a+4] not in pattern:
                        per_itemlist.append(item[a+2:a+6])
                    else:
                        pass
                    a = item.find('DVAD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase2.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase2.append(['XXXX'])
        except IndexError:
            caspase2.append(['XXXX'])

    #Caspase 3
    caspase3 = []
    per_itemlist = []
    pattern = ['P','E','D','Q','K','R']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('DMQD',1)
            if x_count == 0:
                caspase3.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('DMQD', a, -1)           
                    if item[a+4] not in pattern:
                        per_itemlist.append(item[a+2:a+6])
                    else:
                        pass
                    a = item.find('DMQD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase3.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase3.append(['XXXX'])
        except IndexError:
            caspase3.append(['XXXX'])

    #Caspase 4
    caspase4 = []
    per_itemlist = []
    pattern = ['P','E','D','Q','K','R']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('LEVD',1)
            if x_count == 0:
                caspase4.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('LEVD', a, -1)           
                    if item[a+4] not in pattern:
                        per_itemlist.append(item[a+2:a+6])
                    else:
                        pass
                    a = item.find('LEVD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase4.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase4.append(['XXXX'])
        except IndexError:
            caspase4.append(['XXXX'])

    #Caspase 5

    caspase5 = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            x_count = item.count('EHD',1)
            if x_count == 0:
                caspase5.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('EHD', a, -1)           
                    if (item[a-1] == 'L') or (item[a-1] == 'W'):
                        per_itemlist.append(item[a+1:a+5])
                    else:
                        pass
                    a = item.find('EHD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase5.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase5.append(['XXXX'])
        except IndexError:
            caspase5.append(['XXXX'])
    
    # Caspase 6
    caspase6 = []
    per_itemlist = []
    pattern = ['P','E','D','Q','K','R']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('VE',1)
            if x_count == 0:
                caspase6.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('VE', a, -1)           
                    if (item[a+2] == 'H') or (item[a+2] == 'I'):
                        if item[a+3] == 'D':
                            if item[item+4] not in pattern:
                                per_itemlist.append(item[a+2:+6])
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                    a = item.find('VE',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase6.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase6.append(['XXXX'])
        except IndexError:
            caspase6.append(['XXXX'])

    #Caspase 7
    caspase7 = []
    per_itemlist = []
    pattern = ['P','E','D','Q','K','R']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('DEVD',1)
            if x_count == 0:
                caspase7.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('DEVD', a, -1)           
                    if item[a+4] not in pattern:
                        per_itemlist.append(item[a+2:+6])
                    else:
                        pass
                    a = item.find('DEVD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase7.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase7.append(['XXXX'])
        except IndexError:
            caspase7.append(['XXXX'])
    
    #Caspase8
    caspase8 = []
    per_itemlist = []
    pattern = ['P','E','D','Q','K','R']
    for item in PeptideLibraryList:
        try:
            x_count = item.count('ETD',1)
            if x_count == 0:
                caspase8.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('EHD', a, -1)           
                    if (item[a-1] == 'I') or (item[a-1] == 'L'):
                        if item[a+3] not in pattern:
                            per_itemlist.append(item[a+1:a+5])
                        else:
                            pass
                    else:
                        pass
                    a = item.find('EHD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase8.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase8.append(['XXXX'])
        except IndexError:
            caspase8.append(['XXXX'])
        
    #Caspase 9
    caspase9 = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            x_count = item.count('LEHD',1)
            if x_count == 0:
                caspase9.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('LEHD', a, -1)   
                    per_itemlist.append(item[a+2:a+6])        
                    a = item.find('EHD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase9.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase9.append(['XXXX'])
        except IndexError:
            caspase9.append(['XXXX'])

    #Caspase 10
    caspase10 = []
    per_itemlist = []
    for item in PeptideLibraryList:
        try:
            x_count = item.count('IEAD',1)
            if x_count == 0:
                caspase10.append(['XXXX'])
            else:
                a = 1
                while x_count > 0:   
                    a = item.find('IEAD', a, -1)   
                    per_itemlist.append(item[a+2:a+6])        
                    a = item.find('IEAD',a+1,-1)
                    x_count = x_count - 1

                if len(per_itemlist) > 0:
                    caspase10.append(per_itemlist)
                    per_itemlist = []
                else:
                    caspase10.append(['XXXX'])
        except IndexError:
            caspase10.append(['XXXX'])

    Misc_Peptidases_Zip = []
    Misc_Peptidases_Zip = [ list(e) for e in zip(Arg_C, Asp_N, BNPS, Clost, CNBr, Formic_Acid, Glutamyl, Iodoso, LysC, NTCB, Enterokinase, FactorXa, GranzymeB, Elastase, Hydroxylamine, Proline_Endopeptidase, ProteinaseK, Thermolysin, Staph, Thrombin, caspase1, caspase2, caspase3, caspase4, caspase5, caspase6, caspase7, caspase8, caspase9, caspase10)]
    
    ## print(len(Arg_C), len(Asp_N), len(Asp_N), len(BNPS),len(Clost),len(CNBr),len(Formic_Acid),len(Glutamyl),len(Iodoso),len(NTCB),len(Enterokinase),len(FactorXa),len(GranzymeB),len(Elastase),len(Hydroxylamine),len(Proline_Endopeptidase),len(ProteinaseK),len(Thermolysin),len(Staph),len(Thrombin),len(caspase1),len(caspase2),len(caspase3),len(caspase4),len(caspase5),len(caspase6),len(caspase7),len(caspase8),len(caspase9), len(caspase10))
    
    Misc_Peptidases_Len = []
    templist = []
    for item in Misc_Peptidases_Zip:
        for x in item:
            if x[0] == 'XXXX':
                templist.append(0)
            else:
                templist.append(len(x))
            
        Misc_Peptidases_Len.append(templist)
        templist = []

    Proteases = ["Arg-C:", "Asp-N", "BNPS-Skatole", "Clostripain", "CNBr", "Formic Acid", "Glutamyl", "Iodosobenzoic acid", "LysC", "NTCB", "Enterokinase", "FactorXa", "GranzymeB", "Neutrophil Elastase", "Hydroxylamine", "Proline Endopeptidase", "ProteinaseK", "Thermolysin", "Staphylococcal I", "Thrombin", "Caspase1", "Caspase2", "Caspase3", "Caspase4", "Caspase5", "Caspase6", "Caspase7", "Caspase8", "Caspase9", "Caspase10"]
    Misc_Peptidases_Name_Count = []
  
    for item in Misc_Peptidases_Len:
        temp = zip(Proteases, item)
        Misc_Peptidases_Name_Count.append(list(map(list, temp)))
        

    index_list = []
    for item in Misc_Peptidases_Name_Count:
        for index, value in enumerate(item):
            if value[1] == 0:
                index_list.append(index)
            else:
                continue
           
        for x in sorted(index_list, reverse=True):
            del item[x]
        index_list = []

    return Misc_Peptidases_Name_Count


def OpenReadingFrame(DNASequence,minLen,maxLen):
    '''
    Converts input DNA sequence into a peptide sequence
    Starts translating to peptide, until it runs out of triplets.
    Once translated, returns as tuple: (CDS, peptide)
    '''
    #setting up codon table (could use the translate function from biopython, but want to be able to change this to reflect reprogramming)
    translation = {
    "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
    "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
    "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
    "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",

    "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
    "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
    "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
    "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

    "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
    "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
    "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
    "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",

    "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
    "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
    "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
    "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"}
    
    
    #setting start codon to look for
    #StartCodon = 'ATG'
    
    #making copy of input (need to trim to start codon without affecting final output indexing)
    SubString = DNASequence
    #print 'starting string: ' +str(SubString)

    #looping over every available start codon, from the start, until acceptable peptide found
    #while StartCodon in SubString:
#        print 'found start codon'
        #find first start codon
        #StartIndex = SubString.find(StartCodon)
        #trim sequence to start at start codon
        #SubString = SubString[SubString.find(StartCodon):]
        #print 'current string: ' +str(SubString)
        #set up empty peptide sequence
    PeptideSequence = Initiator
        #making copy of input (need to trim during translation without affecting position of start codon search)
    translationString = SubString[3:]
        #scanning down every three letter set
    while len(translationString) > 2:
#            print 'current string: ' +str(translationString[0:3])+' '+str(translationString[3:])
#            print 'translated to: '+ str(translation[translationString[0:3]])
            #adding codon translation to peptide sequence, or 'N' if there are unclear nucleotides
        if 'N' in translationString[0:3]:
            PeptideSequence += 'X'
        else:
            PeptideSequence += translation[translationString[0:3]]
            #trimming sequence by last translated codon
        translationString = translationString[3:]
#            print 'current peptide: '+str(PeptideSequence)
            #if upper length limit hit, stop translating
            #if (len(PeptideSequence)>maxLen):
#                print '--length limit hit--'
            #    return SubString[0:(len(PeptideSequence)*3)],PeptideSequence
            #    break
            #if in frame stop codon found, check if length is ok. Either accept and stop or reject and continue
            #if '+' in PeptideSequence:
            #    if len(PeptideSequence)>minLen:
            #        return SubString[0:(len(PeptideSequence)*3)],PeptideSequence
            #        break
            #    else:
            #        break
            #elif '*' in PeptideSequence:
            #    if len(PeptideSequence)>minLen:
            #        return SubString[0:(len(PeptideSequence)*3)],PeptideSequence
            #        break
            #    else:
            #        break
        #remove first 3 bases(start codon) and continue search
        #SubString = SubString[3:]
    #print PeptideSequence
    if ('*' in PeptideSequence) or ('+' in PeptideSequence):
        return (DNASequence,PeptideSequence)
        #only generates an output sequence if there is an in fram stop codon
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# return a list of lists with peptide-sequences and their frequencies, sorted by frequency in descending order
def SortedPeptideSequencesList(fastqFileLocation,minLen,maxLen):
    '''
    Opens a single sequence file and processes those into a sorted abundance list
    Finds exact matches to forward and reverse primers, and passses the sequence between these to the translation function
    This means the translation function doesn't need to find the right start codon
    '''
    #opening the file, reading out the text, closing the file
    RawDataFile = open(fastqFileLocation, 'r')
    Lines = RawDataFile.readlines()
    RawDataFile.close

    #defining empty list and dictionary
    PeptideSequences = {}
    PeptideSequencesList = []
    
    #using the last character before the file extension to define the round number
    SelectionRoundNumber = fastqFileLocation[fastqFileLocation.find('.')-1]
    
    #taking input file name without .fastq extension and appending text and new extension type
    CSVFileName1 = fastqFileLocation[:-6] + '_allreads_rnd'+str(SelectionRoundNumber)+'.csv'
    ResultsSummaryFile1 = open(CSVFileName1, 'w')
    
    #setting up output file with column headers
    ResultsSummaryFile1.write('ORF' + ',')
    ResultsSummaryFile1.write('peptide')
    ResultsSummaryFile1.write('\n')
    
    for Line in Lines:
        #only taking through sequences with exact matches to forward and reverse primers (may exclude some sequences with poor sequence quality or single misreads)
        if ('ATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATG' in Line) and ('TAGGACGGGGGGCGGGAGGCGGG' in Line):
            startIndex = Line.find('ATACTAATACGACTCACTATAGGATTAAGGAGGTGATATTTATG') + 45
            endIndex = Line.find('TAGGACGGGGGGCGGGAGGCGGG') + 3
            #trimming to only coding sequence
            Line = Line[startIndex:endIndex]
            #print Line
            #translating resulting sequence, returns a cds and peptide sequence as a tuple
            Line = OpenReadingFrame(Line,minLen,maxLen)
            if (Line != None) and (len(Line[1])>minLen) and (len(Line[1])<maxLen):
                #writing tuple resulting from translation to output file
                ORF = Line[0]
                ResultsSummaryFile1.write(str(ORF)+',')
                #print str(ORF)
                Peptide = Line[1]
                ResultsSummaryFile1.write(str(Peptide)+'\n')
                #print str(Peptide)
                #incrementing peptide in count file (dictionary), or creating new entry if not yet present
                if Peptide not in PeptideSequences:
                    PeptideSequences[str(Peptide)] = 1
                else:
                    PeptideSequences[str(Peptide)] = PeptideSequences[str(Peptide)] + 1
                   
    ResultsSummaryFile1.close()
    
    # convert the dictionary into a list of lists
    for key, value in PeptideSequences.items():
        PeptideSequencesList.append([str(SelectionRoundNumber), key, value])
    # sort the PeptideSequenceList by peptide sequence occurence in descending order
    SortedPeptideSequences = sorted(PeptideSequencesList, key = lambda x: x[2], reverse = True)
    return SortedPeptideSequences
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

def SelectionResultsSummary(DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName):
    '''
    Opens multiple files and organises the data resulting from each into the final overview files
    '''
    import os #allows interaction with operating system
    #minLen = input("Minimum sequence length desired: ")
    #maxLen = input("Maximum sequence length desired: ")
    minLen = 18
    maxLen = 100
    
    # (1) create ConcatenatedResultsList using Results from all the rounds of selection       
    ConcatenatedResultsList = {}
    # open files in directory one by one
    for file in os.listdir(DataFolderLocation):
        CurrentFile = ''
        if file.endswith('.fastq'): # this conditional is necessary; without it some shit appears in the beginning of the file list
            CurrentFile = os.path.join(DataFolderLocation, file)
    # print CurrentFile
            
    # (1.A) extract round number from the file name        
            SelectionRoundNumber = CurrentFile[CurrentFile.find('.')-1]
    # print SelectionRoundNumber
    
    # (1.B) extract single round results        
            SingleRoundResults = SortedPeptideSequencesList(CurrentFile,minLen,maxLen)
                    
    # (1.C) populate ConcatenatedResultsList                
            ConcatenatedResultsList[SelectionRoundNumber] = SingleRoundResults
    # print ConcatenatedResultsList
    
    # create sorted list of rounds
    RoundsList = sorted(ConcatenatedResultsList.keys())
    # print RoundsList
    
    # create list of total number of sequences in round
    TotalNumberOfPeptidesList = []
    for Round in RoundsList:
        RoundData = ConcatenatedResultsList[Round]
        NumberOfPeptidesInARound = 0
        for PeptideInstance in range(len(RoundData)):
            NumberOfPeptidesInARound = NumberOfPeptidesInARound + RoundData[PeptideInstance][2]
        TotalNumberOfPeptidesList = TotalNumberOfPeptidesList + [NumberOfPeptidesInARound]
    print("sequences extracted per round: "+ str(TotalNumberOfPeptidesList))
    
    # create list of TopNPeptidesNumber (N) from BaseSelectionRoundNumber (K)
    TopNPeptidesList = []
    for Number in range(TopNPeptidesNumber):
        TopNPeptidesList = TopNPeptidesList + [ConcatenatedResultsList[str(BaseSelectionRoundNumber)][Number][1]]
    # print TopNPeptidesList
    
    # create an empty list of lists of occurences of peptides (k) in every round (n) - for peptide absolute numbers
    ListOfPeptidesOccurancesByRound = []
    # create an empty list of lists of fractions of peptides (k) in every round (n) - for peptide percentages
    ListOfPeptidesFractionsByRound = []
    
    for SelectedPeptide in TopNPeptidesList:
    # create an empty list of peptide occurences by round (n)
        PeptideOccurancesByRound = []
    # create an empty list of peptide fractions by round (n)
        PeptideFractionsByRound = []
        for Round in RoundsList:
            AbundanceInARound = 0
            FractionInARound = 0
            RoundData = ConcatenatedResultsList[Round]
            for PeptideData in RoundData:
                PeptideSequence = PeptideData[1]
                if PeptideSequence == SelectedPeptide:
                    AbundanceInARound = PeptideData[2]
                    FractionInARound = float(PeptideData[2])/float(TotalNumberOfPeptidesList[RoundsList.index(Round)])
            PeptideOccurancesByRound = PeptideOccurancesByRound + [AbundanceInARound]
            PeptideFractionsByRound = PeptideFractionsByRound + [FractionInARound]
        ListOfPeptidesOccurancesByRound = ListOfPeptidesOccurancesByRound + [PeptideOccurancesByRound]
    # print ListOfPeptidesOccurancesByRound
        ListOfPeptidesFractionsByRound = ListOfPeptidesFractionsByRound + [PeptideFractionsByRound]
    # print ListOfPeptidesFractionsByRound
    #-------------------------------------------------------------------------------
    
    #-------------------------------------------------------------------------------
    #Generating lists of peptide parameters
    #truncating all sequences by 7 to remove an assumed GSGSGS* linker (not needed for parameters, but keep in sequence list to find truncations)
    if expanded == 1:
        peptideParameters = []
        for item in TopNPeptidesList:
            peptideParameters.append(str(item[:-7]))

        ChargepH7 = CheckNetCharge(7,peptideParameters)
        Looplist = Loops_Size(peptideParameters)
        LipophilicNumber = AmountLip(peptideParameters)
        hydrophobs = GRAVY(peptideParameters)
        Extinc280 = E_280(peptideParameters)
        Extinc214 = E_214(peptideParameters)
        cellperm = CellPermiability(peptideParameters)
        isoelec = pI(peptideParameters)
        TrypsinCuts = Trypsin(TopNPeptidesList)
        PepsinLow = Pepsin_pH1_3(TopNPeptidesList)
        PepsinHigh = Pepsin_pH3(TopNPeptidesList)
        ChymoCut = Chymotrypsin_LowSpecificity(TopNPeptidesList)
        MiscCuts = Misc_peptidases(TopNPeptidesList)

    #generating output file
    #-------------------------------------------------------------------------------
    CSVFileName = SummaryFileName + '.csv'
    SelectionResultsSummaryFile = open(CSVFileName, 'w')
    
    SelectionResultsSummaryFile.write('peptide sequence' + ',')
    for Round in RoundsList:
        SelectionResultsSummaryFile.write('round # ' + Round + ' occurence (#)' + ',')
    if expanded == 1:
        SelectionResultsSummaryFile.write('Sequence length' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Molecular weight (av)' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Charge at pH 7' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Number of loops' +','+ 'Size of first loop'+','+ 'Size of second loop' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Fraction of lipophilic residues' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('GRAVY hydrophobicity' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('E280' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('E214' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Predicted cell permeability' +',' '(hydrophobicity, polarity)' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Charge transition pH (positive to neutral)' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Charge transition pH (neutral to negative)' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Trypsin cut sites' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Pepsin cut sites (pH 1.3)' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Pepsin cut sites (pH 3)' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Chymotrypsin cut sites' + ',') #adds column header (duplicate and change string to add another)
        SelectionResultsSummaryFile.write('Number of other protease cut sites' + ',') #adds column header (duplicate and change string to add another)
    SelectionResultsSummaryFile.write('\n')
    
    for PeptideIndex in range(len(TopNPeptidesList)):
        SelectionResultsSummaryFile.write(TopNPeptidesList[PeptideIndex] + ',')
        SummaryRoundsData = ListOfPeptidesOccurancesByRound[PeptideIndex]
        for RoundDatum in SummaryRoundsData:
            SelectionResultsSummaryFile.write(str(RoundDatum) + ',')
        if expanded == 1:
            SelectionResultsSummaryFile.write(str(peptideLength(TopNPeptidesList[PeptideIndex])) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(MolecularWeight(peptideParameters[PeptideIndex])) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(ChargepH7[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(Looplist[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(LipophilicNumber[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(hydrophobs[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(Extinc280[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(Extinc214[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(cellperm[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(isoelec[PeptideIndex][0]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(isoelec[PeptideIndex][1]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(TrypsinCuts[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(PepsinLow[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(PepsinHigh[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(ChymoCut[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
            SelectionResultsSummaryFile.write(str(MiscCuts[PeptideIndex]) + ',') #adds output of new function (duplicate and change function call to add another)
        SelectionResultsSummaryFile.write('\n')
        
    SelectionResultsSummaryFile.write('total #' + ',')
    for DataByRound in TotalNumberOfPeptidesList:
        SelectionResultsSummaryFile.write(str(DataByRound) + ',')
    SelectionResultsSummaryFile.write('\n\n\n')
    
    SelectionResultsSummaryFile.write('peptide sequence' + ',')
    for Round in RoundsList:
        SelectionResultsSummaryFile.write('round # ' + Round + ' fraction (%)' + ',')
    SelectionResultsSummaryFile.write('\n')
    
    for PeptideIndex in range(len(TopNPeptidesList)):
        SelectionResultsSummaryFile.write(TopNPeptidesList[PeptideIndex] + ',')
        SummaryRoundsData = ListOfPeptidesFractionsByRound[PeptideIndex]
        for RoundDatum in SummaryRoundsData:
            SelectionResultsSummaryFile.write('{:.3%}'.format(RoundDatum) + ',')
        SelectionResultsSummaryFile.write('\n')
            
    SelectionResultsSummaryFile.close()
    #-------------------------------------------------------------------------------
    #creating a plot of abundance changes across rounds
    #-------------------------------------------------------------------------------
    #import matplotlib.pyplot as plt
    
    #for PeptideNumber in range(len(TopNPeptidesList)):
    #    plt.plot(RoundsList, ListOfPeptidesFractionsByRound[PeptideNumber])
    
    #plt.xlabel('Selection Round #', fontsize=14)
    #plt.ylabel('Peptide Fraction', fontsize=14)
    #legend = plt.legend(TopNPeptidesList, loc='upper center', bbox_to_anchor=(0.5, -0.15),
    #        fancybox=True, shadow=False, ncol=2)
    #PNGFileName = SummaryFileName + '.png'
    #plt.savefig(PNGFileName, bbox_extra_artists=[legend],bbox_inches='tight', dpi = 300)
    #plt.show()
#-------------------------------------------------------------------------------

#_____________________________RUNNING THE FUNCTION_____________________________#

#change the file name and directory below to the desired output file name and the directory where the unarchived fastq files are 
#(make sure last character of the file name is the round number)
#change the 'base selection round number' to the selection round (from file name) that the peptides will be ranked on (usually the last)
#change the 'top peptides number' to the number of peptides to be plotted
#___DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName___
#SelectionResultsSummary('/Users/seino/surfdrive/python/Visual Studio Code testing/Sample sequences/', 4, 100, '/Users/seino/surfdrive/python/Visual Studio Code testing/Sample sequences/test')
SelectionResultsSummary(DataFolderLocation, BaseSelectionRoundNumber, TopNPeptidesNumber, SummaryFileName)
print('Done')
