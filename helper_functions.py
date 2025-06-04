def insertIntoString(string,char,position):
    return string[:position] + char + string[position:len(string)]

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('----the brown cat', 'thes--e brownies-', 5.0)

    Other alignments are also possible.

    """
    lengthSeq1 = len(seq1)
    lengthSeq2 = len(seq2)

    optMatrix = [[0 for j in range(lengthSeq2 + 1)] for i in range(lengthSeq1 + 1)]
    tracebackMatrix = [[('ins') if i == 0 else ('del') for j in range(lengthSeq2 + 1)] for i in range(lengthSeq1 + 1)]
    tracebackMatrix[0][0] = ('')
    for i in range(1, lengthSeq1 + 1):
        optMatrix[i][0] = i * scoring_function(seq1[i - 1],"-") 
    for j in range(1, lengthSeq2 + 1):
        optMatrix[0][j] = j * scoring_function("-",seq2[j - 1]) 

    for i in range(1, lengthSeq1 + 1):
        for j in range(1, lengthSeq2 + 1):
            allign = optMatrix[i-1][j-1] + scoring_function(seq1[i -1],seq2[j-1])
            delete = optMatrix[i - 1][j] + scoring_function("-", seq2[j - 1])
            insert = optMatrix[i][j-1] + scoring_function(seq1[i - 1],"-")
            optMatrix[i][j] = max(allign, delete, insert)
            
            if optMatrix[i][j] == allign:
                tracebackMatrix[i][j] = ('align')  
            elif optMatrix[i][j] == insert:
                tracebackMatrix[i][j] = ('ins') 
            elif optMatrix[i][j] == delete:
                tracebackMatrix[i][j] = ('del')    


    tracebackSteps = []
    nextOperation = tracebackMatrix[lengthSeq1][lengthSeq2]
    xCounter = lengthSeq1
    yCounter = lengthSeq2
    while True:
        tracebackSteps.append(nextOperation)
        if nextOperation== 'align':
            xCounter -= 1
            yCounter -= 1
        elif nextOperation == 'ins':
            yCounter -= 1
        elif nextOperation == 'del':
            xCounter -= 1
        nextOperation = tracebackMatrix[xCounter][yCounter]
        if nextOperation == (''):
            break
    tracebackSteps = tracebackSteps[::-1]

    for i in range(0,len(tracebackSteps)):
        action = tracebackSteps[i]
        if action == "allign":
            pass
        elif action == "ins":
            seq1 = insertIntoString(seq1,"-",i)
        elif action == "del":
            seq2 = insertIntoString(seq2,"-",i)

    return (seq1, seq2, optMatrix[lengthSeq1][lengthSeq2]) 


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('the-- brown', 'these brown', 7.0)

    Other alignments are also possible.

    """
    lengthSeq1 = len(seq1)
    lengthSeq2 = len(seq2)

    optMatrix = [[0 for j in range(lengthSeq2 + 1)] for i in range(lengthSeq1 + 1)]
    
    tracebackMatrix = [[('ins') if i == 0 else ('del') for j in range(lengthSeq2 + 1)] for i in range(lengthSeq1 + 1)]
    tracebackMatrix[0][0] = ('')
    for i in range(1, lengthSeq1 + 1):
        optMatrix[i][0] = max(i * scoring_function(seq1[i - 1],"-"),0) 
    for j in range(1, lengthSeq2 + 1):
        optMatrix[0][j] = max(j * scoring_function("-",seq2[j - 1]),0) 

    maxScore = -1
    maxScoreIndexS1 = 0
    maxScoreIndexS2 = 0
    for i in range(1, lengthSeq1 + 1):
        for j in range(1, lengthSeq2 + 1):
            allign = optMatrix[i-1][j-1] + scoring_function(seq1[i -1],seq2[j-1])
            delete = optMatrix[i - 1][j] + scoring_function("-", seq2[j - 1])
            insert = optMatrix[i][j-1] + scoring_function(seq1[i - 1],"-")
            elt = optMatrix[i][j] = max(allign, delete, insert, 0)
            
            if elt == allign or 0:
                tracebackMatrix[i][j] = ('align')  
            elif elt == insert:
                tracebackMatrix[i][j] = ('ins') 
            elif elt == delete:
                tracebackMatrix[i][j] = ('del')                 
            if elt > maxScore:
                maxScore = elt
                maxScoreIndexS1 = i
                maxScoreIndexS2 = j

    tracebackSteps = []
    nextOperation = tracebackMatrix[lengthSeq1][lengthSeq2]
    xCounter = maxScoreIndexS1
    yCounter = maxScoreIndexS2
    while True:
        tracebackSteps.append(nextOperation)
        if nextOperation == 'align':
            xCounter -= 1
            yCounter -= 1
        elif nextOperation == 'ins':
            yCounter -= 1
        elif nextOperation == 'del':
            xCounter -= 1
        nextOperation = tracebackMatrix[xCounter][yCounter]
        if nextOperation == ('') or optMatrix[xCounter][yCounter] == 0:
            break
    tracebackSteps = tracebackSteps[::-1]

    seq1 = seq1[xCounter:maxScoreIndexS1]
    seq2 = seq2[yCounter:maxScoreIndexS2]
    for i in range(0,len(tracebackSteps)):
        action = tracebackSteps[i]
        if action == "allign":
            pass
        elif action == "ins":
            seq1 = insertIntoString(seq1,"-",i)
        elif action == "del":
            seq2 = insertIntoString(seq2,"-",i)


    return (seq1, seq2, optMatrix[maxScoreIndexS1][maxScoreIndexS2])
