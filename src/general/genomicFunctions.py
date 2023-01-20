def convert_IUPAC_to_regular_expression(nucleotides):
    iupacToRegExDict = {
        "A":"A",
        "C":"C",
        "G":"G",
        "T":"T",
        "W":"[AT]",
        "S":"[CG]",
        "R":"[AG]",
        "Y":"[CT]",
        "K":"[GT]",
        "M":"[AC]",
        "B":"[CGT]",
        "D":"[AGT]",
        "H":"[ACT]",
        "V":"[ACG]",
        "N":"[ACGT]"
    }
    regEx = ""
    for n in nucleotides:
        regEx += iupacToRegExDict[n]
    return regEx
