import sys, getopt

def Get_argv(argv):
    try:
        options, args = getopt.getopt(argv, 'hi:o:')
    except getopt.GetoptError:
        print ('Error while reading arguments. FormatDESeqInput.py -h')
        sys.exit(2)
    if len(options) == 0:
        print ('Error while reading arguments. FormatDESeqInput.py -h')
        sys.exit(2)
    for option, arg in options:
        if len(options) == 1:
            if option == '-h':
                print ('Merge/format expression level counts from multiple samples.')
                print ('Usage: FormatDESeqInput.py -i <input file1>:<input file2>:<input file3>:... -o <output file>')
                sys.exit()
            else:
                print ('Error while reading arguments. FormatDESeqInput.py -h')
                sys.exit(2)
        elif len(options) == 2:
            if option == '-i':
                if arg == '':
                    print ('Error while reading arguments. FormatDESeqInput.py -h')
                    sys.exit(2)
                input_file = arg.split(':')
                print('Input files: '+', '.join(input_file))
            elif option == '-o':
                if arg == '':
                    print ('Error while reading arguments. FormatDESeqInput.py -h')
                    sys.exit(2)
                output_file = arg
                print('Output file: '+output_file)
                
    return input_file, output_file


def Read_coverage(infile):
    try:
        temp_dic = {}
        with open(infile) as temp:
            for i in temp:
                j = i.split('\t')
                temp_dic[j[8]] = j[9]

    except FileNotFoundError:
        print ('No such file or directory: '+infile)
        sys.exit(2)

    return temp_dic                


def Write_coverage(outfile, coverages, order):
    CDSs = list(coverages[order[0]].keys())
    CDSs.sort()
    try:
        with open(outfile, 'w') as temp:
            temp.write('\t'.join(['GI', ]+[i.split('/')[-1][:-4] for i in order]))
            temp.write('\n')
            for CDS in CDSs:
                temp.write(CDS)
                for sample in order:
                    temp.write('\t'+coverages[sample][CDS])
                temp.write('\n')
                
    except FileNotFoundError:
        print ('No such file or directory: '+outfile)
        sys.exit(2)

    return False


args = Get_argv(sys.argv[1:])

coverages = {}
for sample in args[0]:
    coverages[sample] = Read_coverage(sample)

Write_coverage(args[1], coverages, args[0])
