# File: FASTA_converter_1.13.py


## This function, true to it's name, performs the main processes of the program. It reads the BED input file, while writing the output FASTA file.

def Main_process(BEDin_in, FASTAin_in, RnaDNA_in, OUTPUTin):     
        BEDin = BEDin_in
        FASTAin = FASTAin_in
        OUTPUTraw3 = RnaDNA_in


        BED = open(BEDin, "r")
        OUTPUT = open(OUTPUTin, "w")

        line = BED.readline()

        words = line.split()


        if words[0][0] == '#':
                 line = BED.readline()
        else:
                A = 0


        previous_word = 'NA'
        Match = 0
        Match1st = 0
        SignMinus = "-"

        while line:
            words = line.split()
                           
            print( ">" + words[0] + " - " + words[1] + " - " + words[2] + " - " + words[3])
            OUTPUT.write( ">" + words[0] + " - " + words[1] + " - " + words[2] + " - " + words[3] + '\n')
            current_word = words[0]
            start = int(words[1])
            end = int(words[2])
            more_info = words[3]
            signout = words[4]
    

            if previous_word == current_word:
                 output_print( line2[start:end], signout, SignMinus, OUTPUTraw3,OUTPUT)
                          
            else:
                FASTA = open(FASTAin, "r")
                line2 = FASTA.readline()
    
                while line2:
                       if words[0] in  line2 :
                            line2 = FASTA.readline()
                            words2 = line2.split()
                
                            output_print( line2[start:end], signout, SignMinus,OUTPUTraw3,OUTPUT)
                          
                            FASTA.close()
                            Match = 1
                            break
                       else:
                            line2 = FASTA.readline()
                            words2 = line2.split()
            
               

            previous_word = current_word    
            line = BED.readline()

        print(' ')


        print("Output saved.")
        print(OUTPUTin)
        print(' ')

        print( "© Sohan Choudhury")
        print( "Department of Hematology | Yale University")

        BED.close()
        OUTPUT.close()

        print(' ')



## This function does two things. If the BED input line has a "-", it will find the complementary nucleotides of the corresponding sequence, then reverse it.
        ## If the user asks for the output FASTA file to be RNA, the function will replace all the T's in the corresponding sequence with U's.

def output_print(string_in, signout_in, SignMinus_in, OUTPUTraw3_in, OUTPUT):   
        if signout_in == SignMinus_in:
            dict = str.maketrans("ATGC", "TACG")
            value = string_in
            result = value.translate(dict)
            if OUTPUTraw3_in == "y":
                 print( result[::-1].replace("T","U"))
                 OUTPUT.write(result[::-1].replace("T","U") + '\n')
            else:
                 print( result[::-1])
                 OUTPUT.write(result[::-1] + '\n')
        else:
            if OUTPUTraw3_in == "y":
                 print( string_in.replace("T","U"))
                 OUTPUT.write( string_in.replace("T","U") + '\n')
            else:            
                print( string_in)
                OUTPUT.write(string_in + '\n')



## Below is the user interface, created with argparse. The first two arguments, the input BED and FASTA files, are required. The rest are optional.

import os, argparse

parser = argparse.ArgumentParser(description='This is a bioinformatics file format converter that utilizes Python 3. Input your BED file as well as a FASTA sequence to obtain the respective genome sequence in the FASTA format. Order does not matter. Only first two fields are required. Input format: programlocation -b BEDfilelocation -f Sequencelocation -r y/n -o Output')
parser.add_argument("-b", required = True, type=argparse.FileType('r'),   help="BED file *required*")
parser.add_argument("-g", required = True, type=argparse.FileType('r'),  help="Genome sequence file *required*")
parser.add_argument("-r", type=str, help="RNA for output FASTA? (y/n)" , choices=['y', 'n'], nargs = '?', default = 'n')
parser.add_argument("-o", type=argparse.FileType('w'), help="Name of output file", nargs = '?')

parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
args = parser.parse_args()

if args.verbose:
    print ("Arguments:  ".format(args.g,  args.b, args.r, args.o))

else:

    if args.o:
            OutputFileName = args.o.name
    else:
            
            OutputFileName = args.b.name.replace("bed", "fa")
            

    Main_process( args.b.name, args.g.name, args.r, OutputFileName)
