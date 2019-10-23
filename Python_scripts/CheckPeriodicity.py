import time

t1=time.time()

import sys, getopt

def Get_argv(argv):
    try:
        options, args = getopt.getopt(argv, 'hi:a:o:')
    except getopt.GetoptError:
        print ('Error while reading arguments. Call CheckPeriodicity.py -h')
        sys.exit(2)
    if len(options) == 0:
        print ('Error while reading arguments. Call CheckPeriodicity.py -h')
        sys.exit(2)
    for option, arg in options:
        if len(options) == 1:
            if option == '-h':
                print ('Check RiboSeq integrity by periodicity.')
                print ('Usage: CheckPeriodicity.py -i <input BED file> -a <annotation GFF file> -o <output file>')
                sys.exit()
            else:
                print ('Error while reading arguments. Call CheckPeriodicity.py -h')
                sys.exit(2)
        elif len(options) == 3:
            if option == '-i':
                if arg == '':
                    print ('Error while reading arguments. Call CheckPeriodicity.py -h')
                    sys.exit(2)
                input_file = arg
                print('Input BED: '+input_file)
            elif option == '-a':
                if arg == '':
                    print ('Error while reading arguments. Call CheckPeriodicity.py -h')
                    sys.exit(2)
                annot = arg
                print('Input annotation: '+annot)
            elif option == '-o':
                if arg == '':
                    print ('Error while reading arguments. Call CheckPeriodicity.py -h')
                    sys.exit(2)
                output_file = arg
                print('Output file: '+output_file)
    return input_file, annot, output_file


def Parse_annotation(infile):
    try:
        with open(infile) as temp:
            temp2 = []
            for i in temp:
                j = i.split('\n')[0].split('\t')
                temp2.append([j[3], j[4], j[6]])
            
    except FileNotFoundError:
        print ('No such file or directory: '+infile)
        sys.exit(2)
    return temp2


def Generate_profile(infile):
    profile5, profile3 = {}, {}
    try:
        with open(infile) as temp:
            for i in temp:
                j = i.split('\n')[0].split('\t')
                if j[5] == '+':
                    try:
                        profile5[j[5]+str(int(j[1])+1)] += 1
                    except:
                        profile5[j[5]+str(int(j[1])+1)] = 1
                    try:
                        profile3[j[5]+j[2]] += 1
                    except:
                        profile3[j[5]+j[2]] = 1
                else:
                    try:
                        profile5[j[5]+j[2]] += 1
                    except:
                        profile5[j[5]+j[2]] = 1
                    try:
                        profile3[j[5]+str(int(j[1])+1)] += 1
                    except:
                        profile3[j[5]+str(int(j[1])+1)] = 1
    except FileNotFoundError:
        print ('No such file or directory: '+infile)
        sys.exit(2)
    return profile5, profile3

def Calculate_density(profile, CDS):
    if CDS[2] == '+':
        ranges = range(int(CDS[0])-100, int(CDS[0])+301), range(int(CDS[1])-300,int(CDS[1])+101)
    else:
        ranges = range(int(CDS[1])+100,int(CDS[1])-301,-1), range(int(CDS[0])+300, int(CDS[0])-101,-1)
    
    temp = []
    count, num_reads = -100, 0
    for pos in ranges[0]:
        try:
            temp.append(profile[CDS[2]+str(pos)])
            num_reads += profile[CDS[2]+str(pos)]
        except:
            temp.append(0)
        count += 1
    if num_reads < 100:
        return None
    else:
        temp2 = [round(float(i)/float(max(temp)),4) for i in temp]
        return temp2

def Average_density(end_5p, end_3p):
    master_5p, master_3p = [], []
    for pos in range(0,401,1):
        pos_sum = 0
        for CDS in end_5p:
            pos_sum += CDS[pos]
        master_5p.append(round(pos_sum/len(end_5p),4))
        pos_sum = 0
        for CDS in end_3p:
            pos_sum += CDS[pos]
        master_3p.append(round(pos_sum/len(end_3p),4))
    return master_5p, master_3p
    
infile, annot, outfile = Get_argv(sys.argv[1:])

annotation = Parse_annotation(annot)
profiles = Generate_profile(infile)

end_5p, end_3p = [], []
for CDS in annotation:
    if int(CDS[1])-int(CDS[0]) > 300:
        temp = Calculate_density(profiles[0], CDS)
        if temp != None:
            end_5p.append(temp)
        temp2 = Calculate_density(profiles[1], CDS)
        if temp2 != None:
            end_3p.append(temp2)
            
del profiles

avg_5p, avg_3p = Average_density(end_5p, end_3p)

try:
    output = open(outfile, 'w')
    output.write('Position_relative_to_start_codon\t')
    output.write("5'_end_assigned (n=%s)\t3'_end_assigned (n=%s)\n" %(len(end_5p), len(end_3p)))
    for i in range(401):
        output.write('%s\t%s\t%s\n' %(i-100,avg_5p[i],avg_3p[i]))
    output.close()
except FileNotFoundError:
    print ('No such file or directory: '+outfile)
    sys.exit(2)
    
print('% sec' %round(time.time()-t1,4))