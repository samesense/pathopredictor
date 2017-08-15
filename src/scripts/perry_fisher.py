"""Two sided fisher test."""
import sys, myfisher, fisher, pandas, numpy, math, itertools

def main(in_file, out_file):
    df = pandas.read_csv(in_file, sep='\t')
    df['pval'] = myfisher.fisherTestVec(numpy.array(df['pos_fam']),
                                        numpy.array(df['fg_other']),
                                        numpy.array(df['ac']),
                                        numpy.array(df['bg_other']))
    df['fg_gtr'] = (numpy.array(df['pos_fam']) /
                    numpy.array(df['fg_tot'])) > (numpy.array(df['ac']) /
                                                  numpy.array(df['bg_tot']))

    df.to_csv(out_file, index=False, sep='\t')

if __name__ == "__main__":
    in_file, out_file = sys.argv[1:]
    main(in_file, out_file)
