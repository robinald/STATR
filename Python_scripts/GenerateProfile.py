import sys, getopt

def Get_argv(argv):
    try:
        options, args = getopt.getopt(argv, 'hi:e:a:n:o:')
    except getopt.GetoptError:
        print ('Error while reading arguments. GenerateProfile.py -h')
        sys.exit(2)
    if len(options) == 0:
        print ('Error while reading arguments. GenerateProfile.py -h')
        sys.exit(2)
    for option, arg in options:
        if len(options) == 1:
            if option == '-h':
                print ('Generate genome-wide RiboSeq profile by user-selected read end.')
                print ('Usage: GenerateProfile.py -i <input BED file> -e <end, either 5 or 3> -n <normalization factor> -o <output file>')
                sys.exit()
            else:
                print ('Error while reading arguments. GenerateProfile.py -h')
                sys.exit(2)
        elif len(options) == 5:
            if option == '-i':
                if arg == '':
                    print ('Error while reading arguments. GenerateProfile.py -h')
                    sys.exit(2)
                input_file = arg
                print('Input BED: '+input_file)
            elif option == '-e':
                if arg == '':
                    print ('Error while reading arguments. GenerateProfile.py -h')
                    sys.exit(2)
                print("%s' end of the reads will be used" %arg)
                end = arg
            elif option == '-a':
                if arg == '':
                    print ('Error while reading arguments. GenerateProfile.py -h')
                    sys.exit(2)
                acc = arg
            elif option == '-n':
                if arg == '':
                    print ('Error while reading arguments. GenerateProfile.py -h')
                    sys.exit(2)
                norm = arg
            elif option == '-o':
                if arg == '':
                    print ('Error while reading arguments. GenerateProfile.py -h')
                    sys.exit(2)
                output_file = arg
                print('Output file: '+output_file)
                
    return input_file, end, acc, norm, output_file

def Generate_profile(infile, end):
    profile = {}
    read_num = 0
    if end == '5':
        try:
            with open(infile) as temp:
                for i in temp:
                    read_num += 1
                    j = i.split('\n')[0].split('\t')
                    if j[5] == '+':
                        try:
                            profile[j[5]+str(int(j[1])+1)] += 1
                        except:
                            profile[j[5]+str(int(j[1])+1)] = 1
                    else:
                        try:
                            profile[j[5]+j[2]] += 1
                        except:
                            profile[j[5]+j[2]] = 1

        except FileNotFoundError:
            print ('No such file or directory: '+infile)
            sys.exit(2)
            
    elif end == '3':
        try:
            with open(infile) as temp:
                for i in temp:
                    read_num += 1
                    j = i.split('\n')[0].split('\t')
                    if j[5] == '+':
                        try:
                            profile[j[5]+j[2]] += 1
                        except:
                            profile[j[5]+j[2]] = 1
                    else:
                        try:
                            profile[j[5]+str(int(j[1])+1)] += 1
                        except:
                            profile[j[5]+str(int(j[1])+1)] = 1
        except FileNotFoundError:
            print ('No such file or directory: '+infile)
            sys.exit(2)
            
    else:
        print ('Wrong read end. Select either 5 or 3')
        sys.exit(2)
    return profile, float(read_num)


def Get_accession(infile):
    try:
        input_file = open(infile)
        for i in range(10):
            temp = input_file.readline()
        accession = temp.split('\t')[0]
        input_file.close()

    except FileNotFoundError:
        print ('No such file or directory: '+infile)
        sys.exit(2)
    return accession


def Print_profile(profile, sample, acc, outfile, norm, read_num):
    temp_dic = {'+':sample+'_F', '-':sample+'_R'}
    prefix = [acc, 'STATR']
    suffix = ['.', '.\n']
    if norm == 'N':
        try:
            with open(outfile, 'w') as temp:
                for pos, score in profile.items():
                    strand = pos[0]
                    line = prefix + [temp_dic[strand], pos[1:], pos[1:], strand+str(score), strand] +suffix
                    temp.write('\t'.join(line))
        except:
            print ('No such file or directory: '+outfile)
            sys.exit(2)
    else:
        try:
            norm2=float(norm)
    
            nfactor = read_num/norm2
            try:
                with open(outfile, 'w') as temp:
                    for pos, score in profile.items():
                        strand = pos[0]
                        line = prefix + [temp_dic[strand], pos[1:], pos[1:], strand+str(round(float(score)/nfactor,2)), strand] +suffix
                        temp.write('\t'.join(line))
            except:
                print ('No such file or directory: '+outfile)
                sys.exit(2)
    
        except:
            print ('Wrong normalization paramter.')
            sys.exit(2)

    return True


args = Get_argv(sys.argv[1:])
profile, read_num = Generate_profile(args[0], args[1])
Print_profile(profile, args[0][:-4], Get_accession(args[2]), args[4], args[3], read_num)
del profile

