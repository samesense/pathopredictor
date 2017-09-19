"""Parse blat to find genomcic coords for transcript coords."""
import argparse

def main(args):
    with open(args.output, 'w') as fout, open(args.blat_input) as f:
        print('c\tchrom\tg\tc_nuc\tstrand', file=fout)
        for x in range(5):
            f.readline()
        line = f.readline().strip().split()
        strand = line[8]
        chrom = line[13]
        transcript_seq = line[21]
        # (transcript_start_ls, transcript_len_ls,
        #  genome_start_ls) = [ x.strip(',').split(',') for x in line[18:21] ]
        for t_len, t_st, g_st in zip(*[ x.strip(',').split(',') for x in line[18:21] ] ):
            #print(t_st, g_st, t_len)
            for offset in range( int(t_len) ):
                t = int(t_st)+1+offset
                g = int(g_st)+1+offset
#                print( len(transcript_seq), t_st, t_len)
                c_nuc = transcript_seq[int(t_st) + offset].upper()
                print('\t'.join([str(t), chrom, str(g), c_nuc, strand]), file=fout)

if __name__ == "__main__":
    desc = 'Mk coordinate map.'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('blat_input', 'output',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)
