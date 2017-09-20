"""Parse blat to find genomcic coords for transcript coords."""
import argparse

def load_seq(fa_file):
    with open(fa_file) as f:
        f.readline()
        return f.readline().strip()

def write_hash(transcript_seq, line, output_file):
    strand = line[8]
    chrom = line[13]
    #transcript_seq = line[21].replace(',', '')
        # (transcript_start_ls, transcript_len_ls,
        #  genome_start_ls) = [ x.strip(',').split(',') for x in line[18:21] ]
#    print(line[11:13])
    with open(output_file, 'w') as fout:
        print('c\tchrom\tg\tc_nuc\tstrand', file=fout)
        for t_len, t_st, g_st in zip(*[ x.strip(',').split(',') for x in line[18:21] ] ):
                #print(t_st, g_st, t_len)
#            print(t_len, t_st, len(transcript_seq))
            for offset in range( int(t_len) ):
                t = int(t_st)+1+offset
                g = int(g_st)+1+offset
    #                print( len(transcript_seq), t_st, t_len)
                c_nuc = transcript_seq[int(t_st) + offset].upper()
                print('\t'.join([str(t), chrom, str(g), c_nuc, strand]), file=fout)

def main(args):
    transcript_seq = load_seq(args.fa_file)
    
    blat_data = ''
    blat_max = 0
    with open(args.blat_input) as f:
        for x in range(5):
            f.readline()
        for l in f:
            line = l.strip().split()
            query_st, query_end = [int(x) for x in line[11:13]]
            q_size = abs(query_end - query_st)
            if q_size > blat_max:
                blat_max = q_size
                blat_data = line

    # process query with biggest span
    line = blat_data
    write_hash(transcript_seq, line, args.output)

if __name__ == "__main__":
    desc = 'Mk coordinate map.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('blat_input', 'fa_file', 'output',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
