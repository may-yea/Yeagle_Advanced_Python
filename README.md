# Yeagle_Advanced_Python


# Sequence Objects

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

































# Sequence Annotation

```python
# Import the sequence record

from Bio.SeqRecord import SeqRecord
```


```python
# import sequence from Bio.Seq

from Bio.Seq import Seq
```


```python
# define the sequence

simple_seq = Seq("GATC")
```


```python
# create sequence record

simple_seq_r = SeqRecord(simple_seq)
```


```python
# bring the sequence record forward

simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# pass the ID to the record

simple_seq_r.id = "AC12345"
```


```python
# create a description for the sequence

simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
# print the description

print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
# call the sequence

simple_seq_r.seq
```




    Seq('GATC')




```python
# bring the record forward

simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
# store annotations in your record

simple_seq_r.annotations["evidence"] = "None. This is an example"
```


```python
# print annotation

print(simple_seq_r.annotations["evidence"])
```

    None. This is an example



```python
# bring the record forward

simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
# create letter annotations

simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
# print letter annotations

print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# save complete sequence from this link

# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
# import SeqIO function from biopython

from Bio import SeqIO
```


```python
# tell the function to read the saved sequence and tell the function the file type

record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
# pull the record of the saved sequence

record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# pull the sequence

record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# pull the ID

record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
# pull the description

record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# you can only get the information that is maintained in the FASTA file
# part 1: no annotations

record.annotations
```




    {}




```python
# you can only get the information that is maintained in the FASTA file
# part 2: no references

record.dbxrefs
```




    []




```python
# you can only get the information that is maintained in the FASTA file
# part 3: no features

record.features
```




    []




```python
# download and save this file

# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
# tell the SeqIO function to read the genbank file

record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
# open the file

record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# pull the sequence

record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# pull the ID

record.id
```




    'NC_005816.1'




```python
# pull the name

record.name
```




    'NC_005816'




```python
# pull the description

record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# genbank does not provide letter annotations (absent in file)

record.letter_annotations
```




    {}




```python
# pull the number of annotations

len(record.annotations)
```




    13




```python
# pull the source -> gives the species

record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
# pull the reference

record.dbxrefs
```




    ['Project:58037']




```python
# pull the number of features

len(record.features)
```




    41




```python
# import the SeqFeature function from Biopython

from Bio import SeqFeature
```


```python
# tell the function where you think the start position is roughly (after X letter)

start_pos = SeqFeature.AfterPosition(5)
```


```python
# tell the function where you think the sequence ends (ambiguously)

end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
```


```python
# create a location

my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
# print the location

print(my_location)
```

    [>5:(8^9)]



```python
# pull out the start location only

my_location.start
```




    AfterPosition(5)




```python
# pull out the end location only

my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
# pull out just the integer for the end location

int(my_location.end)
```




    9




```python
# pull out only the start location integer

int(my_location.start)
```




    5




```python
# pass the numbers to the verbs

exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
# print the exact location

print(exact_location)
```

    [5:9]



```python
# pull the exact start location

exact_location.start
```




    ExactPosition(5)




```python
# Import empty seqrecord

from Bio.SeqRecord import SeqRecord
```


```python
# create a seqrecord

record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
                       "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
                       "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
                       "SSAC"),
                   id= "gi|14150838|gb|AAK54648.1|AF376133_1",
                   description="chalcone synthase [Cucumis sativus]",
                  )
```


```python
# print the file and clarify the format

print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
# print the record

print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
# Read the sequence

record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
# pull the record

record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# pull the number of base pairs in the record

len(record)
```




    9609




```python
# pull the number of features in the record

len(record.features)
```




    41




```python
# print the X (20th) feature from the record

print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
# print the Xth (21st) feature from the record

print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# sub divide the record within a specified region of the location

sub_record = record[4300:4800]
```


```python
# get the length of the sub record

len(sub_record)
```




    500




```python
# get the length of the sub record features

len(sub_record.features)
```




    2




```python
# look at the first of the two features in the sub record

sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
# look at the second of the two features in the sub record

sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
# print the first feature of the sub record

print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
# print the second feature of the sub record

print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# look at the subrecord annotations

sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
# there are no subrecord references

sub_record.dbxrefs
```




    []




```python
# add to the annotations of the record

sub_record.annotations["topology"] = "linear"
```


```python
# pull the new subrecord annotations

sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
# pull the subrecord ID

sub_record.id
```




    'NC_005816.1'




```python
# pull the subrecord name

sub_record.name
```




    'NC_005816'




```python
# pull the subrecord description

sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# edit the subrecord description

sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
# pull the new description

sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# print everyting through 200 from the subrecord and clarify the format

print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
# go back to using the genbank record

record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
# pull the record

record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# pull the length of the record

len(record)
```




    9609




```python
# pull the number of features of the record

len(record.features)
```




    41




```python
# pull the reference of the record

record.dbxrefs
```




    ['Project:58037']




```python
# pull the annotation keys (titles)

record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# shift the starting position of the circle organism to the end

shifted = record[2000:] + record[:2000]
```


```python
# pull the shifted record

shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# pull the length of the shifted record

len(shifted)
```




    9609




```python
# pull the number of features of the shifted

len(shifted.features)
```




    40




```python
# pull the shifted annotation keys

shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
# pull the shifted reference

shifted.dbxrefs
```




    []




```python
# define the shifted reference the same as the record

shifted.dbxrefs = record.dbxrefs[:]
```


```python
# pull the redefined shifted reference

shifted.dbxrefs
```




    ['Project:58037']




```python
# define the shifted annotations as the record annotations

shifted.annotations = record.annotations.copy()
```


```python
# pull the redefined shifted annotations

shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# pull the record

record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# "print the first value as a string and the next four as integers" => define the values of the record

print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# define the reverse complement

rc = record.reverse_complement(id ="Testing")
```


```python
# pull the reverse complement

rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# print the same definitions of values, but for the reverse complement

print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



```python

```
