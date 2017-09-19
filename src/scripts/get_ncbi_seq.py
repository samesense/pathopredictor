"""Get ncbi seq for transcript"""
import eutils.client, argparse

def main(args):
    ec = eutils.client.Client()
    egs = ec.efetch(db='nuccore', id=args.transcript)
    seq = egs.gbseqs[0].sequence
    with open(args.output, 'w') as fout:
        print('>' + args.transcript, file=fout)
        print(seq, file=fout)

if __name__ == "__main__":
    desc = 'Pull ncbi seq to use in blat'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('transcript', 'output',)
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)

