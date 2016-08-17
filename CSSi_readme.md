# Yale
##Sohan Choudhury

Cryptic Splice Site Identifier. Given an annotation file (either GTF or GFF3) as well as a STAR .tab file containing splice junctions, the program will find all non-canonical splice sites. Program works in all cases: whether 3', 5', or both splice sites are cryptic. Output is printed and saved, contains key information for analysis.

Notable features:
- Multiple .tab files can be used simultaneously, as user has option to input .txt file containing paths of an unlimited number of .tab files
- Output contains both cryptic splice junction as well as nearest annotated case
- Output also provides distance of novel case from canonical case, as well as directionality of each case
- Saved output is tab delimited for hassle-free analysis using spreadsheets such as Excel, also contains number of uniquely mapped reads for each case for easy statistical significance calculations for each novel case
- Parser style user input


