# Yeagle_Advanced_Python

```python
# Load Biopython package

from Bio.Seq import Seq
```


```python
# Input a sequence

my_seq = Seq("GATCG")
```


```python
# Number the sequence

for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# Print the length of the sequence

print(len(my_seq))
```

    5



```python
# Print the letter from a specific position of the sequence

print(my_seq[0])
```

    G



```python
# Print the letter from a specific position of the sequence Pt.2

print(my_seq[4])
```

    G



```python
# Print the letter from a specific position of the sequence Pt.3

print(my_seq[2])
```

    T



```python
# Create a new sequence and count number of occurences of a partial sequence or letter
# Excludes overlapping partials/letters

Seq("AAAA").count("AA")
```




    2




```python
# Create a new sequence

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# Get the length

len(my_seq)
```




    32




```python
# Count the occurences of a letter

my_seq.count("G")
```




    9




```python
# Measure the GC content (% form)

100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    46.875




```python
# Import the Biopython function for determining GC content

from Bio.SeqUtils import gc_fraction
```


```python
# Use the function for determining GC content (decimal form)

gc_fraction(my_seq)
```




    0.46875




```python
# Create a slice of the sequence from the 4th letter to the 12th

my_seq[4:12]
```




    Seq('GATGGGCC')




```python
# Create a slice of a start and a stop in a stride
# Starts at the first base pair and cuts out every third

my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
# Create a slice that starts at the base pair in position 1 and cuts out every third

my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
# We can choose where we initiate a slice and where it ends

my_seq[2:3]
```




    Seq('T')




```python
# Print the entire sequence backward by starting the slice at the opposite end

my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
# Turn a seq object into a string

str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
# Create a placeholder using python string formatting
# We can use this to give sequences labels

fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
# Print the placeholder 

print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
# Create two sequences

seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
# Add sequences together

seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
# We can also add sequences in different orders

seq2 + seq1
```




    Seq('AACCGGACGT')




```python
# Create contigs sequences

contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
# Create spacer with specified length

spacer = Seq("N" *10)
```


```python
# Join the contigs together with the spacer

spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
# Identify a sequence with mixed case

dna_seq = Seq("acgtACGT")
```


```python
# Print the sequence 

dna_seq
```




    Seq('acgtACGT')




```python
# Make the sequence lowercase

dna_seq.lower()
```




    Seq('acgtacgt')




```python
# Make the sequence uppercase

dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# Make the sequence equal to the uppercase sequence

dna_seq = dna_seq.upper()
```


```python
# We cannot find a fragment that is in the wrong case

"gtac" in dna_seq
```




    False




```python
# We can find a fragment by using the proper case

"GTAC" in dna_seq
```




    True




```python
# Create a new sequence

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# Print the complement

my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
# Print the reverse complement

my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
# Create a protein sequence and print the complement

protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
# Create a coding dna sequence

coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
# Print the sequence

coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# Create a template dna sequence
# The template dna should equal to the reverse complement of the coding dna

template_dna = coding_dna.reverse_complement()
```


```python
# Print the template dna sequence

template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
# Reprint the coding dna sequence to compare to the template

coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# Create a messenger rna sequence 
# The messenger rna should equal the transcribed coding dna

messenger_rna = coding_dna.transcribe()
```


```python
# Print the messenger rna

messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# The transribed reverse complement of the template dna should be equal to the messenger rna

template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# We can back transcribe the messenger rna to get the coding dna

messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# Recall messenger rna

messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# Translate messenger rna into protein sequence
# The asterisks indicate stop codons

messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# When we specify the codon table we eliminate the premature stop codon

coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
# We can also use the NCBI table number

coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
# We can also make the protein sequence stop at the initial stop codon

coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
# We can combine the NCBI table number function and initial stop codon function
# to get the full protein sequence without the asterisk

coding_dna.translate(table =2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
# We can modify the stop codon to a different symbol

coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
# Input bacterial gene base pair sequence

gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
# Translate the gene with the bacterial table to get the nucleotide sequence

gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
# Publish all the way to the stop codon without the asterisk

gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# Tell the Biopython that this is a complete coding dna sequence
# so the translation starts with the proper codon (methionine)

gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
# Import a codon table package

from Bio.Data import CodonTable
```


```python
# Use codon table package function (or NBCI table #)
# Input standard table (Table 1)

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
# Use codon table package function (or NBCI table #) Pt.2
# Input mitochondrial table (Table 2)

mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# Print standard codon table

print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# Print vertebrate mitochondrial codon table

print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# We can identify the stop codons

mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
# We can identify the start codons

mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# Create seq_object

seq = Seq("ACGT")
```


```python
# We ask does this sequence equal the sequence we created?

"ACGT" == seq1
```




    True




```python
# Does the sequence we created equal this sequence?
# True to both, they are fully equal

seq1 == "ACGT"
```




    True




```python
# We creat a seq object with unkown letters but a known length

unknown_seq = Seq(None, 10)
```


```python
# When we access the unknown seq object we see there is no sequence
# data other than the length

unknown_seq
```




    Seq(None, length=10)




```python
# We are able to get the length

len(unknown_seq)
```




    10




```python
# Pull the sequence from the file. Identify the starting position of the sequence 
# and the sequence length.

seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# When we try to access the sequence within an undefined range
# we see the sequence data is empty other than the length.

seq[1000:1020]
```




    Seq(None, length=20)




```python
# We can access the sequence data for a range that is defined
# Recall the starting position of our sequence (117512683)

seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# We can have partial information on a sequence combined with unkown information

seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
# Define part of a sequence

seq = Seq("ACGT")
```


```python
# Identify an undefined sequence with a defined length

undefined_seq = Seq(None, length =10)
```


```python
# Combine the partial sequence with the undefined part

seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
# Define a sequence (immutable by default)

my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
# Import mutable sequence package

from Bio.Seq import MutableSeq
```


```python
# Redefine sequence as mutable

mutable_seq = MutableSeq(my_seq)
```


```python
# Print new sequence

mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# We can now edit the sequence by changing position

mutable_seq[5] = "C"
```


```python
# Print the edited sequence

mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# We can also edit the sequence by removing the first of any letter

mutable_seq.remove("T")
```


```python
# Print the edit

mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# We can reverse the sequence

mutable_seq.reverse()
```


```python
# Print

mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# We can redefine our sequence as the new sequence to make the new version
# immutable again

new_seq = Seq(mutable_seq)
```


```python
# Print

new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# If we do not want to se the seq object functions directly we can still perform
# them by importing them for use with strings

from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
# Define string

my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
# Perform reverse complent with string instead of seq object

reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
# Transcribe the string instead of seq object

transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
# Back transcribe the string instead of the seq object

back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
# Translate into protein without using my seq

translate(my_string)
```




    'AVMGRWKGGRAAG*'




```python

```
