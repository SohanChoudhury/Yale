# File: Cryptic3_1.14.py

def File_process(g_file_in):
    
        ArrayElement= []
        
        GFILE = open(g_file_in, "r")
        #OUTPUT = open(OUTPUTin, "w")

        #print(OUTPUTin)
        
        line = GFILE.readline()
        i = 0
        while line:
            
            words = line.split()
            if words[0][0] == '#':
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
            

                ArrayElement[i][5] = gene_id
                ArrayElement[i][6] = transcript_id
                ArrayElement[i][7] = gene_name
                ArrayElement[i][8] = gene_type
                ArrayElement[i][9] = exon_number

                i = i+1
 
            
                line = GFILE.readline()
        
        GFILE.close()
        #OUTPUT.close()
        #print (" Array : " + str(ArrayElement))
        return ArrayElement





def find(g3_array_in,len,chromosome_in,start_in,end_in,signout_in, unique_reads_in, OUTPUT):
    # Try to do exact Match
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
                    #print ( "Chromosome, direction, and start of splice junction matched : " + str(row) + "  " + str(g3_array_in[row]))
                    StartFound = 1
                    exon_start_Next_number = g3_array_in[row][9]
                    exon_start_Next_number_int = int(exon_start_Next_number) +1
                    exon_start_Next_number = str(exon_start_Next_number_int)
                    row = row +1
        else:
            row = row +1


    row = 0
    while row <= len  and EndFound == 0 :
            
            if (g3_array_in[row][0] == chromosome_in) and (g3_array_in[row][2] == end_in ) and (g3_array_in[row][4] == signout_in) and ((g3_array_in[row][9] == exon_start_Next_number) or ( StartFound == 0)) :
                #print ( "Exact end of splice junction found : " +str(row) + "  "+ str(g3_array_in[row]))
                EndFound = 1
                exon_end_prev_number = g3_array_in[row][9]
                exon_end_prev_number_int = int(exon_end_prev_number) - 1
                exon_end_prev_number = str(exon_end_prev_number_int)
                row = row +1
            else:
                row = row +1

    if ( StartFound ==1 ) and ( EndFound == 1):
        #print( "-------- Exact start and end of splice junction found --------")
        b = 0
    elif ( StartFound == 1 ) and ( EndFound == 0):
        #print ( "********* Exact start found, but exact end of splice junction not found *********")                    
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
        PrevABSdistance = 999999
        row_distance = 0
        

        while row <= len-1:
            g3_array_end_in = int(g3_array_in[row][3])
            g3_array_start_in = int(g3_array_in[row][2])
            if (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) and (g3_array_in[row][9] == exon_start_Next_number):
                #print ( "read the same sequence: " +str(row) + "  " + str(g3_array_in[row]))
                row=row +1
                distance =  end_in_new - g3_array_start_in
                
                ABSdistance = abs(distance)
                #print ( "Distance " + str(ABSdistance))
                row_distance = row -1
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

        #print( " end " + str(end_in_new) + " new end " +  g3_array_in[max_row][2] )        
        if (end_in_new - int(int(g3_array_in[max_row][2]))) < 0 and (g3_array_in[max_row][4] == '+'):
            #Direction  = '-'
            Stream = 'upstream'

        elif (end_in_new - int(int(g3_array_in[max_row][2]))) < 0 and (g3_array_in[max_row][4] == '-'):
            #Direction = '-'
            Stream = 'downstream'

        elif (end_in_new - int(int(g3_array_in[max_row][2]))) > 0 and (g3_array_in[max_row][4] == '+'):
            #Direction = '+'
            Stream = 'downstream'
        elif (end_in_new - int(int(g3_array_in[max_row][2]))) > 0 and (g3_array_in[max_row][4] == '-'):
            #Direction = '+'
            Stream = 'upstream'
        else:
            #print( "Failed.")
            #OUTPUT.write( ("Failed.") + '\n')
            y = 10
        #print ( " Direction " +     Direction)
        #print ( " Lowest Distance : " + str(ABSdistance) + " row :" + str(row_distance) + " Direction: " + Direction +" Diff : "+ str(end_in_new - int(int(g3_array_in[row_distance][2]))))       
        #print ( " max end : " + str(max_end) + " end : " + str(end_in_new));        
        if max_end > end_in_new:
           print(g3_array_in[max_row][7] + "; " + g3_array_in[max_row][5] + "    " + g3_array_in[max_row][4] + "    " + chromosome_in + " : " + start_in + " - " + end_in + "    " + g3_array_in[max_row][0] + " : " + start_in + " - " + g3_array_in[max_row][2] + "     " + str(ABSdistance) +
                    "    " + Stream + "    " + unique_reads_in)
           OUTPUT.write(g3_array_in[max_row][7] + "; " + g3_array_in[max_row][5] + "\t" + g3_array_in[max_row][4] + "\t" + chromosome_in + " : " + start_in + " - " + end_in + "\t" + g3_array_in[max_row][0] + " : " + start_in + " - " + g3_array_in[max_row][2] + "\t" + str(ABSdistance) +
                    "\t" + Stream + "\t" + unique_reads_in + '\n')
           #print ( " End is good, max end:" +str(max_end)+ str(row) + "  "+ str(g3_array_in[max_row]))
        else:
           #print ( " End is not good")
            e = 0
          
    elif ( StartFound == 0 ) and ( EndFound == 1):
        #print ( "********* Exact end found, but exact start of splice junction not found *********")                    
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
                row_distance = row -1
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

        

        if (start_in_new - int(int(g3_array_in[min_row][3]))) < 0 and (g3_array_in[min_row][4] == '+'):
            #Direction  = '-'
            Stream = 'upstream'

        elif (start_in_new - int(int(g3_array_in[min_row][3]))) < 0 and (g3_array_in[min_row][4] == '-'):
            #Direction = '-'
            Stream = 'downstream'

        elif (start_in_new - int(int(g3_array_in[min_row][3]))) > 0 and (g3_array_in[min_row][4] == '+'):
            #Direction = '+'
            Stream = 'downstream'
        elif (start_in_new - int(int(g3_array_in[min_row][3]))) > 0 and (g3_array_in[min_row][4] == '-'):
            #Direction = '+'
            Stream = 'upstream'
        else:
            #print( "Failed.")
            OUTPUT.write("Failed." + '\n')
            e = 12


                            
        #print ( " Lowest Distance : " + str(ABSdistance) + " row :" + str(row_distance) + " Direction: " + Direction + " Diff : " + str(start_in_new - int(int(g3_array_in[row_distance][3]))))  
        if min_start < start_in_new:
           print (  g3_array_in[min_row][7] + "; " + g3_array_in[min_row][5] + "    " + g3_array_in[min_row][4] + "    " + chromosome_in + " : " + start_in + " - " + end_in + "    " + g3_array_in[min_row][0] + " : " + g3_array_in[min_row][3] + " - " + end_in + "     " + str(ABSdistance) +
                    "    " + Stream + "    " + unique_reads_in)
           OUTPUT.write(g3_array_in[min_row][7] + "; " + g3_array_in[min_row][5] + "\t" + g3_array_in[min_row][4] + "\t" + chromosome_in + " : " + start_in + " - " + end_in + "\t" + g3_array_in[min_row][0] + " : " + g3_array_in[min_row][3] + " - " + end_in + "\t" + str(ABSdistance) +
                    "\t" + Stream + "\t" + unique_reads_in + '\n')
           #print ( "  Start is good, min start : " +str(min_start)+ str(row) + "  "+ str(g3_array_in[min_row]))
        else:
           #print ( "  Start is not good")
           r = 9

           
      
                
    else:
        # Check Start position
        #print ( " No start or end match found : Last condition")
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
        
        

        while row <= len-1:
            g3_array_start_in = int(g3_array_in[row][2])
            g3_array_end_in = int(g3_array_in[row][3])
            if (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) :
                #print ( " Row:" + str(row) + "  "+ str(g3_array_in[row]))
                row=row +1
                
                distance =  start_in_new - g3_array_end_in
                ABSStartdistance = abs(distance)
                #print ( " Distance " + str(ABSStartdistance))
                row_distance = row -1
                if ABSStartdistance < MinimumStartDistance:
                    MinimumStartDistance = ABSStartdistance
                    row_distance = row -1
                    min_start_row = row -1
                else:
                    ABSStartdistance = ABSStartdistance
                 #print ( " ABS 

                    
                if min_start < g3_array_start_in:
                     min_start = min_start
                else:
                    min_start = g3_array_start_in
                    min_start_row = row -1
 
            else:
                row=row+1
 
        if (start_in_new - int(int(g3_array_in[row_distance][3]))) > 0:
            Direction  = '+'
        else:
            Direction  = '-'



        if (start_in_new - int(int(g3_array_in[min_start_row][3]))) < 0 and (g3_array_in[min_start_row][4] == '+'):
            #Direction  = '-'
            StartStream = 'upstream'

        elif (start_in_new - int(int(g3_array_in[min_start_row][3]))) < 0 and (g3_array_in[min_start_row][4] == '-'):
            #Direction = '-'
            StartStream = 'downstream'

        elif (start_in_new - int(int(g3_array_in[min_start_row][3]))) > 0 and (g3_array_in[min_start_row][4] == '+'):
            #Direction = '+'
            StartStream = 'downstream'
            
        elif (start_in_new - int(int(g3_array_in[min_start_row][3]))) > 0 and (g3_array_in[min_start_row][4] == '-'):
            #Direction = '+'
            StartStream = 'upstream'
        else:
            #print( "Failed.")
            f = 4

        
        #print(" *** exon : " + str(g3_array_in[min_start_row]) + " Min Difference " + str(MinimumStartDistance))                             
        #print ( " *** Lowest Distance : " + str(MinimumStartDistance) + " row :" + str(row_distance) + " Direction: " + Direction + " Diff : " + str(start_in_new - int(int(g3_array_in[row_distance][3]))))  
        
        
        # Check end of spice junction now
        
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
        PrevABSdistance = 999999
        min_end_row = 0

        #print ( " **** Next Exon : " + str(exon_start_Next_number) + "   row " + str(g3_array_in[row]))
        exon_start_Next_number = g3_array_in[min_start_row][9]
        exon_start_Next_number_int = int(exon_start_Next_number) +1
        exon_start_Next_number = str(exon_start_Next_number_int)
        #print ( " **** Next Exon : " + str(exon_start_Next_number))
        

        while row <= len-1:
            g3_array_end_in = int(g3_array_in[row][3])
            g3_array_start_in = int(g3_array_in[row][2])
            if (g3_array_in[row][0] == chromosome_in)  and (g3_array_in[row][4] == signout_in) and (g3_array_in[row][9] == exon_start_Next_number):
                #print ( "read the same sequence: " +str(row) + "  " + str(g3_array_in[row]))
                row=row +1
                distance =  end_in_new - g3_array_start_in
                ABSEnddistance = abs(distance)
                
                #print ( " **** Distance *** " + str(ABSdistance))

                ABSdistance = abs(distance)
                #print ( "Distance " + str(ABSdistance))
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
                
                
        if (end_in_new - int(int(g3_array_in[row_distance][2]))) > 0:
            Direction  = '+'
        else:
            Direction  = '-'

        if (end_in_new - int(int(g3_array_in[max_row][2]))) < 0 and (g3_array_in[max_row][4] == '+'):
            #Direction  = '-'
            EndStream = 'upstream'

        elif (end_in_new - int(int(g3_array_in[max_row][2]))) < 0 and (g3_array_in[max_row][4] == '-'):
            #Direction = '-'
            EndStream = 'downstream'

        elif (end_in_new - int(int(g3_array_in[max_row][2]))) > 0 and (g3_array_in[max_row][4] == '+'):
            #Direction = '+'
            EndStream = 'downstream'
            
        elif (end_in_new - int(int(g3_array_in[max_row][2]))) > 0 and (g3_array_in[max_row][4] == '-'):
            #Direction = '+'
            EndStream = 'upstream'
        else:
            print( "Failed.")
            
        #print ( " *** End *** Lowest Distance : " + str(MinimumDistance) + " row :" + str(row_distance) + " Direction: " + Direction +" Diff : "+ str(end_in_new - int(int(g3_array_in[row_distance][2]))))       
        #print ( "  *** End ***max end : " + str(max_end) + " end : " + str(end_in_new));
        #print ( " *** End row :" + str(g3_array_in[min_row]))

        if min_start < start_in_new:
           #print ( " *** Start is good, min start : " +str(min_start)+ str(row) + "  "+ str(g3_array_in[min_row]))
            print (  g3_array_in[min_start_row][7] + "; " + g3_array_in[min_start_row][5] + "    " + g3_array_in[min_start_row][4] + "    " + chromosome_in + " : " + start_in + " - " + end_in + "    " + g3_array_in[min_start_row][0] + " : " + g3_array_in[min_start_row][3] + " - " + g3_array_in[min_end_row][2] + "     " + str(MinimumStartDistance) + "; " + str(ABSEnddistance) +
                    "    " + StartStream +  "; " + EndStream + "    " + unique_reads_in)
            OUTPUT.write(g3_array_in[min_start_row][7] + "; " + g3_array_in[min_start_row][5] + "\t" + g3_array_in[min_start_row][4] + "\t" + chromosome_in + " : " + start_in + " - " + end_in + "\t" + g3_array_in[min_start_row][0] + " : " + g3_array_in[min_start_row][3] + " - " + g3_array_in[min_end_row][2] + "\t" + str(MinimumStartDistance) + "; " + str(ABSEnddistance) +
                    "\t" + StartStream +  "; " + EndStream + "\t" + unique_reads_in + '\n')
        else:
           #print ( " *** Start is not good")
           t = 2
        



def search_tab( t_file_in, g3_array,len_of_array_in, OUTPUTin):


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
                #print("Direction not found")
                nope = 100

            #print("chromosome : " + chromosome + " start: "+ start + " end : " + end + " sign : " + sign + " Meaning : " + signout)    

            find(g3_array,len_of_array_in,chromosome, start, end, signout, unique_reads, OUTPUT)

            line = TFILE.readline()

        TFILE.close()
        OUTPUT.close()
        

    


import os, argparse

parser = argparse.ArgumentParser(description="This is a cryptic splice site identifier that utilizes Python 3. Input your GFF3/GTF file as well as a STAR .tab output to obtain the correct output. Order of input fields do not matter. Input format: programlocation -g GFF3/GTFfilelocation -t .tabfilelocation")
parser.add_argument("-g", required = True, type=argparse.FileType('r'),   help="GFF3/GTF *required*")
parser.add_argument("-t", required = True, type=argparse.FileType('r'),   help="STAR .tab *required*")
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
        
    GFF3_array = File_process(args.g.name)
    len_of_array =len( GFF3_array)
    search_tab( args.t.name,GFF3_array,len_of_array, OutputFileName)
