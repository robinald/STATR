import sys, getopt
sys.argv

def get_argv(argv):
    input_file = ''
    output_file = ''
    try:
        options, args = getopt.getopt(argv, 'hi:o:')
    except getopt.GetoptError:
        print ('Error while reading argument. Call ParseGenomeAnnotation.py -h')
        sys.exit(2)
    if len(options) == 0:
        print ('Error while reading argument. Call ParseGenomeAnnotation.py -h')
        sys.exit(2)
    for option, arg in options:
        if len(options) == 1:
            if option == '-h':
                print ('ParseGenomeAnnotation.py sweeps unnecessarily long attribute section of GFF3 file downloaded from NCBI.')
                print ('Usage: Parse_genome_annotation.py -i <input_file> -o <output_file>')
                sys.exit()
            else:
                print ('Error while reading argument. Call ParseGenomeAnnotation.py -h')
                sys.exit(2)
        elif len(options) == 2:
            if option == '-i':
                if arg == '':
                    print ('Error while reading argument. Call ParseGenomeAnnotation.py -h')
                    sys.exit(2)
                input_file = arg
                print('Input GFF3: '+input_file)
            elif option == '-o':
                if arg == '':
                    print ('Error while reading argument. Call ParseGenomeAnnotation.py -h')
                    sys.exit(2)
                output_file = arg
                print('Output GFF3: '+output_file)
    return input_file, output_file


def Parse_input(infile):
    try:
        temp = open(infile, 'r')
        temp2 = temp.readline()
        while temp2[0] == '#':
            temp2 = temp.readline()
        temp3 = [temp2]+temp.readlines()
        temp.close()
        temp4 = []
        for i in temp3:
            temp4.append(i.split('\t'))
        accession = temp4[0][0]
        
    except FileNotFoundError:
        print ('No such file or directory: '+infile)
        sys.exit(2)
    return temp4, accession


def Sweep_and_output(input_list, accession, output_file):
    in_count = len(input_list)
    out_count = 0
    output = open(output_file, 'w')
    for i in input_list:
        try:
            if i[2] == 'CDS':
                start, end, strand = i[3], i[4], i[6]
                ID = Get_ID(i[8])
                output.write('\t'.join([accession, 'RefSeq', 'CDS', start, end, '.', strand, '.', ID]) + '\n')
                out_count += 1
        except IndexError:
            pass
    print ('Read %s lines from GFF3' %in_count)
    print ('%s CDSs remained' %out_count)
    output.close()
    return
    
def Get_ID(attribute):
    attributes = attribute.split(';')
    attributes_dic = {}
    for i in attributes:
        attributes_dic[i.split('=')[0]] = i.split('=')[1]
    return attributes_dic['locus_tag']

input_file, output_file =get_argv(sys.argv[1:])
input_list, accession = Parse_input(input_file)
Sweep_and_output(input_list, accession, output_file)