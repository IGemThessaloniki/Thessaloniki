#TRIGGER_1 = hsa-miR-143-3p, TRIGGER_2 = hsa-miR-30e-5p, TRIGGER_3 = hsa-miR-1246, TRIGGER_ 4 = hsa-miR-30e-3p
trigger = ["UGAGAUGAAGCACUGUAGCUC", "UGUAAACAUCCUUGACUGGAAG", "AAUGGAUUUUUGGAGCAGG","CUUUCAGUCGGAUGUUUACAGC"]

#RBS
RBS = "AGAGGAGA"
#linker
linker = "AACCUGGCGGCAGCGCAAAAG"
#GFP
GFP = "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa".upper()
GFP = GFP.replace("T","U")


#write "ViennaRNA" for analysis with ViennaRNA
#"Nupack" for analysis with nupack
selection ="Nupack"

RBS.upper()
linker.upper()
GFP.upper()

paired = 9

unpaired = 9
