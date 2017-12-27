# File: Cryptic3_2.12.py

# This program will identify cryptic 3' and 5' splice sites when given a STAR .tab output (or a text file that contains locations of multiple .tab files) as well as a reference annotation file (GTF or GFF3)
# Created by Sohan Choudhury for Pillai Lab @ Yale University Department of Hematology


def FileMerge(MultipleTAB_in, tabOutFile_out):
# In the case that a text file with multiple .tab files is provided, this funcion is used to merge and sort the .tab files

  import csv
  import shutil

  OUTPUTin = tabOutFile_out
  MultipleTAB = open(MultipleTAB_in,"r")
  
  line = MultipleTAB.readline()
  words = line.split()
  File1 = words[0]
  FirstFile = File1
  File1_old = words[0] +"_old"
  shutil.copy(File1, File1_old)
  FirstRead = 0

  while line:

         if  FirstRead == 1:
                 shutil.copy(OUTPUTin, File1)
                 FirstRead = 1
         line = MultipleTAB.readline()

         if len(line) == 0:
                 output.close()
                 in_file1.close()
                 in_file2.close()
                 shutil.copy(File1_old, FirstFile)
                 break

         output = open(OUTPUTin,"w")
         rowmatch = 0
         firstTime = 0
         FirstRead = 1
         words = line.split()
         File2 =words[0]
         in_file1 = open(File1,"r")
         reader1 = csv.reader((in_file1), delimiter="\t")

         for row1 in reader1:
                y1 = row1[0], row1[1], row1[2], row1[3]
                in_file2 = open(File2,"r")
                reader2 = csv.reader((in_file2), delimiter="\t")

                if firstTime == 1:
                    output.write( "\n")
                output.write(row1[0]+"\t" + row1[1]+"\t" + row1[2]+"\t" + row1[3]+"\t" + row1[4]+ "\t" + row1[5]+"\t" + row1[6])
                firstTime = 1

                for row2 in reader2:
                        z = row2[0], row2[1], row2[2], row2[3]

                        if tuple(z) in [tuple(y1)]:
                                  output.write( "," + row2[6])
                                  rowmatch = 1

         output.write( "\n")
         rowmatch = 1
         in_file1.close()
         in_file2.close()
         in_file1 = open(File2,"r")
         reader1 = csv.reader((in_file1), delimiter="\t")

         for row1 in reader1:
                rowmatch = 0
                y1 = row1[0], row1[1], row1[2], row1[3]
                in_file2 = open(File1,"r")
                reader2 = csv.reader((in_file2), delimiter="\t")

                for row2 in reader2:
                    z = row2[0], row2[1], row2[2], row2[3]

                    if tuple(z) in [tuple(y1)]:
                            rowmatch = 1

                if rowmatch == 0:
                    output.write(row1[0]+"\t" + row1[1]+"\t" + row1[2]+"\t" + row1[3]+"\t" + row1[4]+ "\t" + row1[5]+"\t" + row1[6] + "\n")

         output.close()
         in_file1.close()
         in_file2.close()

         # This is used to sort the merged file
         MergeSort = 'MergeSort.tab'
         outputMergeSort = open(MergeSort,"w")

         with open(OUTPUTin, 'r') as r:

              for line in sorted(r):
                      outputMergeSort.write(line+"")

         outputMergeSort.close()
         shutil.copy(MergeSort, OUTPUTin)


def File_process(g_file_in):
# This function processes the input annotation file (.gtf or .gff3) and puts the neccesary information into an array for easy access 

        import os

        file_name, extension = os.path.splitext(g_file_in)
        ArrayElement= []
        GFILE = open(g_file_in, "r")
        line = GFILE.readline()
        i = 0

        while line:

            words = line.split('\t')

            if words[0][0] == '#' or  words[2] != 'exon':
                 line = GFILE.readline()

            else:
                ArrayElement.append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])
                ArrayElement[i].append([])

                feature_name      = words[0]
                feature_type      = words[2]
                feature_start     = words[3]
                feature_end       = words[4]
                feature_direction = words[6]
                feature_desc      = words[8]
                ArrayElement[i][0] = feature_name
                ArrayElement[i][1] = feature_type
                ArrayElement[i][2] = feature_start
                ArrayElement[i][3] = feature_end
                ArrayElement[i][4] = feature_direction                        
                                 

                if str(extension) == ".gff3":
                # Checks to see if file is GFF3

                  gene_id_pos = feature_desc.find( 'gene_id=')+8
                  gene_id_end = feature_desc.find( ';',gene_id_pos)
                  gene_id = feature_desc[gene_id_pos:gene_id_end]

                  transcript_id_pos = feature_desc.find( 'transcript_id=')+14
                  transcript_id_end = feature_desc.find( ';',transcript_id_pos)
                  transcript_id = feature_desc[transcript_id_pos:transcript_id_end]

                  gene_name_pos = feature_desc.find( 'gene_name=')+10
                  gene_name_end = feature_desc.find( ';',gene_name_pos)
                  gene_name = feature_desc[gene_name_pos:gene_name_end]

                  gene_type_pos = feature_desc.find( 'gene_type=')+10
                  gene_type_end = feature_desc.find( ';',gene_type_pos)
                  gene_type = feature_desc[gene_type_pos:gene_type_end]

                  exon_number_pos = feature_desc.find( 'exon_number=')+12
                  exon_number_end = feature_desc.find( ';',exon_number_pos)
                  exon_number = feature_desc[exon_number_pos:exon_number_end]



                elif str(extension) == ".gtf":
                # Alternatively, checks to see if file is GTF

                  gene_id_pos = feature_desc.find( 'gene_id "')+9
                  gene_id_end = feature_desc.find( '";',gene_id_pos)
                  gene_id = feature_desc[gene_id_pos:gene_id_end]

                  transcript_id_pos = feature_desc.find( 'transcript_id "')+15
                  transcript_id_end = feature_desc.find( '";',transcript_id_pos)
                  transcript_id = feature_desc[transcript_id_pos:transcript_id_end]

                  gene_name_pos = feature_desc.find( 'gene_name "')+11
                  gene_name_end = feature_desc.find( '";',gene_name_pos)
                  gene_name = feature_desc[gene_name_pos:gene_name_end]

                  gene_type_pos = feature_desc.find( 'gene_type "')+11
                  gene_type_end = feature_desc.find( '";',gene_type_pos)
                  gene_type = feature_desc[gene_type_pos:gene_type_end]

                  exon_number_pos = feature_desc.find( 'exon_number ')+12
                  exon_number_end = feature_desc.find( ';',exon_number_pos)
                  exon_number = feature_desc[exon_number_pos:exon_number_end]



                ArrayElement[i][5] = gene_id
                ArrayElement[i][6] = transcript_id
                ArrayElement[i][7] = gene_name
                ArrayElement[i][8] = gene_type
                ArrayElement[i][9] = exon_number

                i = i+1


                line = GFILE.readline()                

        GFILE.close()
        return ArrayElement


def find(g3_array_in,len,chromosome_in,start_in,end_in,signout_in, unique_reads_in, OUTPUT):
# Main processing to identify cryptic splice sites

    # First, tries to do exact match
    row = 0
    StartFound = 0
    EndFound = 0
    len = len -1
    search = 0
    MiddleStartFound = 0
    MiddleEndFound = 0
    exon_start_Next_number  = 0

    while row <= len and StartFound == 0:

        if (g3_array_in[row][0] == chromosome_in) and (g3_array_in[row][3] == start_in ) and (g3_array_in[row][4] == signout_in):
                    StartFound = 1
                    exon_start_Next_number = g3_array_in[row][9]
                    exon_start_Next_number_int = int(exon_start_Next_number) +1
                    exon_start_Next_number = str(exon_start_Next_number_int)
                    start_row_found = row
                    row = row +1

        else:
            row = row +1

    row = 0

    while row <= len  and EndFound == 0 :

            if (g3_array_in[row][0] == chromosome_in) and (g3_array_in[row][2] == end_in ) and (g3_array_in[row][4] == signout_in) and ((g3_array_in[row][9] == exon_start_Next_number) or ( StartFound == 0)) :
                EndFound = 1
                end_row_found = row
                exon_end_prev_number = g3_array_in[row][9]
                exon_end_prev_number_int = int(exon_end_prev_number) - 1
                exon_end_prev_number = str(exon_end_prev_number_int)
                row = row +1

            else:
                row = row +1

    if ( StartFound ==1 ) and ( EndFound == 1):
        a = 0

    elif ( StartFound == 1 ) and ( EndFound == 0):
    # In this case: exact start of splice junction found in annotation file, but exact end not found
        row = 0
        start_in_int = 0
        end_in_int   = 0
        g3_array_end_in = int(g3_array_in[row][3])
        start_in_new    = int(start_in)
        g3_array_start_in = int(g3_array_in[row][2])
        end_in_new    = int(end_in)
        max_end = 0
        distance = 999999999999999
        ABSdistance = 0
        PrevABSdistance = 99999999999
        row_distance = 0

        while row <= len-1:
            g3_array_end_in = int(g3_array_in[row][3])
            g3_array_start_in = int(g3_array_in[row][2])

            if ( g3_array_start_in > start_in_new ) and (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) and (g3_array_in[row][9] == exon_start_Next_number):
                row=row +1
                distance =  end_in_new - g3_array_start_in
                ABSdistance = abs(distance)

                if ABSdistance < PrevABSdistance:
                    ABSdistance = ABSdistance
                    row_distance = row -1

                else:
                    ABSdistance = PrevABSdistance

                PrevABSdistance = ABSdistance    

                if max_end < g3_array_end_in:
                    max_end = g3_array_end_in
                    max_row = row -1

                else:
                    max_end = max_end

            else:
                row=row+1

        if (end_in_new - int(int(g3_array_in[row_distance][2]))) < 0 and (g3_array_in[row_distance][4] == '+'):
            Stream = 'upstream'

        elif (end_in_new - int(int(g3_array_in[row_distance][2]))) < 0 and (g3_array_in[row_distance][4] == '-'):
            Stream = 'downstream'

        elif (end_in_new - int(int(g3_array_in[row_distance][2]))) > 0 and (g3_array_in[row_distance][4] == '+'):
            Stream = 'downstream'

        elif (end_in_new - int(int(g3_array_in[row_distance][2]))) > 0 and (g3_array_in[row_distance][4] == '-'):
            Stream = 'upstream'

        else:
            print( "Failed.")

        if max_end > end_in_new:
           print(g3_array_in[start_row_found][7] + "; " + g3_array_in[start_row_found][5] + "    " + g3_array_in[start_row_found][4] + "    " + chromosome_in + " : " + start_in + " - " + end_in + "    " + g3_array_in[start_row_found][0] + " : " + start_in + " - " + g3_array_in[row_distance][2] + "     " + str(ABSdistance) +
                    "    " + Stream + "    " + unique_reads_in)
           OUTPUT.write(g3_array_in[start_row_found][7] + "; " + g3_array_in[start_row_found][5] + "\t" + g3_array_in[start_row_found][4] + "\t" + chromosome_in + " : " + start_in + " - " + end_in + "\t" + g3_array_in[start_row_found][0] + " : " + start_in + " - " + g3_array_in[row_distance][2] + "\t" + str(ABSdistance) +
                    "\t" + Stream + "\t" + unique_reads_in + '\n')

        else:
            a = 0

    elif ( StartFound == 0 ) and ( EndFound == 1):
    # In this case: exact end of splice junction found in annotation file, but exact start not found    
        row = 0
        start_in_int = 0
        end_in_int   = 0
        g3_array_end_in = int(g3_array_in[row][3])
        start_in_new    = int(start_in)
        g3_array_start_in = int(g3_array_in[row][2])
        start_in_new    = int(start_in)
        end_in_new    = int(end_in)
        min_start = 999999999999999
        min_row = 0
        distance = 999999999999999
        ABSdistance = 0
        row_distance = 0

        while row <= len-1:
            g3_array_start_in = int(g3_array_in[row][2])
            g3_array_end_in = int(g3_array_in[row][3])

            if (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) and (g3_array_in[row][9] == exon_end_prev_number):
                row=row +1
                distance =  start_in_new - g3_array_end_in
                ABSdistance = abs(distance)

                if ABSdistance > g3_array_start_in:
                    ABSdistance = g3_array_start_in
                    row_distance = row -1

                else:
                    ABSdistance = ABSdistance

                if min_start < g3_array_start_in:
                     min_start = min_start

                else:
                    min_start = g3_array_start_in
                    min_row = row -1

            else:
                row=row+1

        if (start_in_new - int(int(g3_array_in[row_distance][3]))) > 0:
            Direction  = '+'

        else:
            Direction  = '-'

        if (start_in_new - int(int(g3_array_in[row_distance][3]))) < 0 and (g3_array_in[row_distance][4] == '+'):
            Stream = 'upstream'

        elif (start_in_new - int(int(g3_array_in[row_distance][3]))) < 0 and (g3_array_in[row_distance][4] == '-'):
            Stream = 'downstream'

        elif (start_in_new - int(int(g3_array_in[row_distance][3]))) > 0 and (g3_array_in[row_distance][4] == '+'):
            Stream = 'downstream'

        elif (start_in_new - int(int(g3_array_in[row_distance][3]))) > 0 and (g3_array_in[row_distance][4] == '-'):
            Stream = 'upstream'

        else:
            print( "Failed.")
                           
        if min_start < start_in_new:
           print (  g3_array_in[row_distance][7] + "; " + g3_array_in[row_distance][5] + "    " + g3_array_in[row_distance][4] + "    " + chromosome_in + " : " + start_in + " - " + end_in + "    " + g3_array_in[row_distance][0] + " : " + g3_array_in[end_row_found][3] + " - " + end_in + "     " + str(ABSdistance) +
                    "    " + Stream + "    " + unique_reads_in)
           OUTPUT.write(g3_array_in[row_distance][7] + "; " + g3_array_in[row_distance][5] + "\t" + g3_array_in[row_distance][4] + "\t" + chromosome_in + " : " + start_in + " - " + end_in + "\t" + g3_array_in[row_distance][0] + " : " + g3_array_in[end_row_found][3] + " - " + end_in + "\t" + str(ABSdistance) +
                    "\t" + Stream + "\t" + unique_reads_in + '\n')

        else:
           print ( "  Start is not good")
               
    else:
    # No start or end match found (both 3' and 5' splice sites are cryptic)
        row = 0
        start_in_int = 0
        end_in_int   = 0
        g3_array_end_in = int(g3_array_in[row][3])
        start_in_new    = int(start_in)
        g3_array_start_in = int(g3_array_in[row][2])
        start_in_new    = int(start_in)
        end_in_new    = int(end_in)
        min_start = 999999999999999
        min_start_row = 0
        distance = 999999999999999
        ABSStartdistance = 0
        row_distance = 0
        MinimumStartDistance = 99999999999999
        goodStartFound = 0
        row_start_distance = 0

        while row <= len-1:
            g3_array_start_in = int(g3_array_in[row][2])
            g3_array_end_in = int(g3_array_in[row][3])

            if (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) :
                row=row +1
                distance =  start_in_new - g3_array_end_in
                ABSStartdistance = abs(distance)

                if ABSStartdistance < MinimumStartDistance:   
                    MinimumStartDistance = ABSStartdistance
                    goodStartFound = 1
                    row_start_distance = row -1

                else:
                    ABSStartdistance = ABSStartdistance

                if min_start < g3_array_start_in:
                     min_start = min_start

                else:
                    min_start = g3_array_start_in
                    min_start_row = row -1

            else:
                row=row+1

        if (start_in_new - int(int(g3_array_in[row_start_distance][3]))) > 0:
            Direction  = '+'

        else:
            Direction  = '-'

        if (start_in_new - int(int(g3_array_in[row_start_distance][3]))) < 0 and (g3_array_in[row_start_distance][4] == '+'):
            StartStream = 'upstream'

        elif (start_in_new - int(int(g3_array_in[row_start_distance][3]))) < 0 and (g3_array_in[row_start_distance][4] == '-'):
            StartStream = 'downstream'

        elif (start_in_new - int(int(g3_array_in[row_start_distance][3]))) > 0 and (g3_array_in[row_start_distance][4] == '+'):
            StartStream = 'downstream'

        elif (start_in_new - int(int(g3_array_in[row_start_distance][3]))) > 0 and (g3_array_in[row_start_distance][4] == '-'):
            StartStream = 'upstream'

        else:
            print( "Failed.")

        row = 0
        start_in_int = 0
        end_in_int   = 0
        g3_array_end_in = int(g3_array_in[row][3])
        start_in_new    = int(start_in)
        g3_array_start_in = int(g3_array_in[row][2])
        end_in_new    = int(end_in)
        max_end = 0
        distance = 999999999999999
        MinimumEndDistance = 999999999999
        ABSEnddistance = 0
        row_distance = 0
        PrevABSdistance = 9999999999
        min_end_row = 0
        max_row = 0

        if goodStartFound == 1:                
                exon_start_Next_number = g3_array_in[row_start_distance][9]
                exon_start_Next_number_int = int(exon_start_Next_number) +1
                exon_start_Next_number = str(exon_start_Next_number_int)

        else:
                a = 0

        while row <= len-1 and goodStartFound == 1:
            g3_array_end_in = int(g3_array_in[row][3])
            g3_array_start_in = int(g3_array_in[row][2])

            if ( g3_array_start_in > start_in_new) and (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) and (g3_array_in[row][9] == exon_start_Next_number):
                row=row +1
                distance =  end_in_new - g3_array_start_in
                ABSEnddistance = abs(distance)
                ABSdistance = abs(distance)
                row_distance = row -1

                if ABSEnddistance < PrevABSdistance:
                    ABSEnddistance = ABSdistance
                    row_distance = row -1
                    min_end_row = row -1

                else:
                    ABSEnddistance = PrevABSdistance

                PrevABSdistance = ABSEnddistance

                if max_end < g3_array_end_in:
                    max_end = g3_array_end_in
                    max_row = row -1

                else:
                    max_end = max_end

            else:
                row=row+1
                
                
        if (end_in_new - int(int(g3_array_in[min_end_row][2]))) > 0:
            Direction  = '+'

        else:
            Direction  = '-'

        if (end_in_new - int(int(g3_array_in[min_end_row][2]))) < 0 and (g3_array_in[min_end_row][4] == '+'):
            EndStream = 'upstream'

        elif (end_in_new - int(int(g3_array_in[min_end_row][2]))) < 0 and (g3_array_in[min_end_row][4] == '-'):
            EndStream = 'downstream'

        elif (end_in_new - int(int(g3_array_in[min_end_row][2]))) > 0 and (g3_array_in[min_end_row][4] == '+'):
            EndStream = 'downstream'
            
        elif (end_in_new - int(int(g3_array_in[min_end_row][2]))) > 0 and (g3_array_in[min_end_row][4] == '-'):
            EndStream = 'upstream'

        else:
            print( "Failed.")
            

        if min_start < start_in_new and max_end > end_in_new:
            print (  g3_array_in[row_start_distance][7] + "; " + g3_array_in[row_start_distance][5] + "    " + g3_array_in[row_start_distance][4] + "    " + chromosome_in + " : " + start_in + " - " + end_in + "    " + g3_array_in[row_start_distance][0] + " : " + g3_array_in[row_start_distance][3] + " - " + g3_array_in[min_end_row][2] + "     " + str(MinimumStartDistance) + "; " + str(ABSEnddistance) +
                    "    " + StartStream +  "; " + EndStream + "    " + unique_reads_in)
            OUTPUT.write(g3_array_in[row_start_distance][7] + "; " + g3_array_in[row_start_distance][5] + "\t" + g3_array_in[row_start_distance][4] + "\t" + chromosome_in + " : " + start_in + " - " + end_in + "\t" + g3_array_in[row_start_distance][0] + " : " + g3_array_in[row_start_distance][3] + " - " + g3_array_in[min_end_row][2] + "\t" + str(MinimumStartDistance) + "; " + str(ABSEnddistance) +
                    "\t" + StartStream +  "; " + EndStream + "\t" + unique_reads_in + '\n')

        else:
            a = 0


def search_tab( t_file_in, g3_array, len_of_array_in, OUTPUTin):
# This function analyzes the .tab file and makes slight adjestments to match the style of the annotation file

        OUTPUT = open(OUTPUTin, "w")
        TFILE = open(t_file_in, "r")
        line = TFILE.readline()
        words = line.split()

        while line:
            words = line.split()
            chromosome = words[0]
            start     = words[1]
            end       = words[2]
            sign      = words[3]
            unique_reads = words[6] 
            startInt = int(start) -1
            endInt   = int(end) -1
            start = str(startInt)
            end   = str(endInt)

            if sign == '2':
                signout = '-'

            elif sign == '1':
                signout = '+'

            else:
                a = 0    

            find(g3_array,len_of_array_in,chromosome, start, end, signout, unique_reads, OUTPUT)

            line = TFILE.readline()

        TFILE.close()
        OUTPUT.close()
        

import csv
import os, argparse

parser = argparse.ArgumentParser(description="This is a cryptic splice site identifier that utilizes Python 3. Input your GFF3/GTF file as well as a STAR .tab output to obtain the correct output. Order of input fields do not matter. Input format: programlocation -g GFF3/GTFfilelocation -t .tabfilelocation")
parser.add_argument("-g", required = True, type=argparse.FileType('r'),   help="GFF3/GTF *required*")
parser.add_argument("-t", required = False, type=argparse.FileType('r'), nargs = '?', help="STAR .tab *required*")
parser.add_argument("-f", required = False, type=argparse.FileType('r'),   help="Text file w/ .tab file locations") 
parser.add_argument("-o", type=argparse.FileType('w'), help="Name of output file", nargs = '?')

parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
args = parser.parse_args()

if args.verbose:
    print ("Arguments:  ".format(args.g))

else:

    if args.o:
        OutputFileName = args.o.name

    else:
        OutputFileName = args.t.name.replace("tab", "tab1")

    if args.f:
    # Checks for multiple .tab files
                 GFF3_array = File_process(args.g.name)
                 len_of_array =len( GFF3_array)
                 tabOutFile_out = 'MergeTab.tab'
                 FileMerge(args.f.name, tabOutFile_out)
                 search_tab( tabOutFile_out,GFF3_array,len_of_array, OutputFileName)

    elif args.t:
                 GFF3_array = File_process(args.g.name)
                 len_of_array =len( GFF3_array)
                 search_tab( args.t.name,GFF3_array,len_of_array, OutputFileName)

    else:
          print ( " One argument( t or f ) is required ")
