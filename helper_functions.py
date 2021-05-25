def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

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
    def scoring_function(x, y):
    if x == "-" or y == "-":
        return -2
    elif x != y:
        return -1
    else:
        return 2
    """

    seqq1="-"+seq1
    seqq2="-"+seq2
    s = [0] * len(seqq1)
    for i in range(len(seqq1)):
        s[i] = [0] * len(seqq2)
    s1 = ""
    s2 = ""
    emp="-"
    for i in range(len(seqq1)):
        for j in range(len(seqq2)):
            if i > 0 and j > 0:
                s[i][j] = max(s[i - 1][j - 1] + scoring_function(seqq1[i], seqq2[j]),
                              s[i - 1][j] + scoring_function(seqq1[i], "-"),
                              s[i][j - 1] + scoring_function("-", seqq2[j]))
            elif i > 0:
                s[i][j] = s[i - 1][j] + scoring_function(seqq1[i], "-")
            elif j > 0:
                s[i][j] = s[i][j - 1] + scoring_function("-", seqq2[j])
            else:
                s[i][j] = 0
    # trackback:
    trackback = []
    i = len(seqq1) - 1
    j = len(seqq2) - 1
    while (i >= 0 and j >= 0):  # za 0 0 da vidam so se desava 2* ka ke e 0 nov while i na soodvetnoto stavam crtki samo
        s[i]
        if i > 0 and j > 0:
            if s[i][j]-scoring_function(seqq1[i], seqq2[j])==s[i-1][j-1]:
                s1 = seqq1[i]+s1
                s2 = seqq2[j]+s2
            elif s[i][j]-scoring_function(seqq1[i], "-")==s[i-1][j]:
                s2 = "-"+s2
                s1 = seqq1[i]+s1
                j=j+1
            elif s[i][j]-scoring_function("-", seqq2[j])==s[i][j-1]:
                s1 = "-"+s1
                s2 = seqq2[j]+s2
                i=i+1
            else:
                print("something is wrong")

        elif i > 0:
            #print("something is wrong i")
            s2 = "-" + s2
            s1 = seqq1[i] + s1
            j = j +1
        elif j > 0:
            #print("something is wrong j")
            s1 = "-" + s1
            s2 = seqq2[j] + s2
            i = i+1
        #else:
            #print("something is wrong 0")
            #s[i][j] = 0

        i = i - 1
        j = j - 1
   
    #s1 = s1[::-1]
    #s2 = s2[::-1]
    #for i in range(len(s)):
     #   print(s[i][:])
    return s1, s2, s[len(seqq1)-1][len(seqq2)-1]



def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

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

    """

    seqq1="-"+seq1
    seqq2="-"+seq2
    s = [0] * len(seqq1)
    for i in range(len(seqq1)):
        s[i] = [0] * len(seqq2)
    s1 = ""
    s2 = ""
    emp="-"
    mx=0
    mi=0
    mj=0
    for i in range(len(seqq1)):
        for j in range(len(seqq2)):
            if i > 0 and j > 0:
                s[i][j] = max(s[i - 1][j - 1] + scoring_function(seqq1[i], seqq2[j]),0,
                              s[i - 1][j] + scoring_function(seqq1[i], "-"),
                              s[i][j - 1] + scoring_function("-", seqq2[j]))
            elif i > 0:
                s[i][j] = max((s[i - 1][j] + scoring_function(seqq1[i], "-")),0)
            elif j > 0:
                s[i][j] = max((s[i][j - 1] + scoring_function("-", seqq2[j])),0)
            else:
                s[i][j] = 0
            if s[i][j]>mx:
                mx=s[i][j]
                mi=i
                mj=j
    # trackback:
    trackback = []
    i = mi#len(seqq1) - 1
    j = mj#len(seqq2) - 1
    while (i >= 0 and j >= 0):  # za 0 0 da vidam so se desava 2* ka ke e 0 nov while i na soodvetnoto stavam crtki samo
        s[i]
        if i > 0 and j > 0:
            if (s[i][j]-scoring_function(seqq1[i], seqq2[j])==s[i-1][j-1]) and s[i][j]>0:
                s1 = seqq1[i]+s1
                s2 = seqq2[j]+s2
            elif (s[i][j]-scoring_function(seqq1[i], "-")==s[i-1][j]) and s[i][j]>0:
                s2 = "-"+s2
                s1 = seqq1[i]+s1
                j=j+1
            elif (s[i][j]-scoring_function("-", seqq2[j])==s[i][j-1]) and s[i][j]>0:
                s1 = "-"+s1
                s2 = seqq2[j]+s2
                i=i+1
            #else:
                #print("probably we have reached the end diagonally :) ")
                #i=0
                #j=0
                

        elif i > 0 and s[i][0]>0:
            #print("something is wrong i")
            s2 = "-" + s2
            s1 = seqq1[i] + s1
            j = j +1
        elif j > 0 and s[0][j]>0:
            #print("something is wrong j")
            s1 = "-" + s1
            s2 = seqq2[j] + s2
            i = i+ 1
        #else:
            #print("probably we have reached the end left/up  may be a mistake :) ")
            #s[i][j] = 0

        i = i - 1
        j = j - 1
   
    #s1 = s1[::-1]
    #s2 = s2[::-1]
    #for i in range(len(s)):
     #   print(s[i][:])
    return s1, s2, mx