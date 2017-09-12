import argparse
import svgutils.transform as sg

def main(args):
    fig = sg.SVGFigure("25cm", "20cm")
    fig1 = sg.fromfile(args.gene_counts)
    fig2 = sg.fromfile(args.mpc_dist)
    fig3 = sg.fromfile(args.fg_roc)
    fig4 = sg.fromfile(args.clinvar_roc)

    plot1 = fig1.getroot()
    plot2 = fig2.getroot()
    plot3 = fig3.getroot()
    plot4 = fig4.getroot()

    plot1.moveto(0, 20, scale=0.75)
    plot2.moveto(400, 0, scale=0.75)
    plot3.moveto(30, 400, scale=0.75)
    plot4.moveto(410, 400, scale=0.75)

    # add text labels
    txt1 = sg.TextElement(25, 20, "A", size=12, weight="bold")
    txt2 = sg.TextElement(400, 20, "B", size=12, weight="bold")

    txt3 = sg.TextElement(25, 390, "C", size=12, weight="bold")
    txt4 = sg.TextElement(400, 390, "D", size=12, weight="bold")

    fig.append([plot1, plot2, plot3, plot4])
    fig.append([txt1, txt2, txt3, txt4])
    fig.save(args.out)
    

    # Figure("10cm", "10cm",
    #        Panel(
    #         Image(120, 120, args.gene_counts),
    #         Text("A", 25, 20, size=12, weigth='bold')
    #         ),
    #        Panel(
    #         Image(120, 120, args.mpc_dist),
    #         Text("B", 25, 20, size=12, weight='bold')
    #         ).move(120, 0),
    #        Panel(
    #         Image(120, 120, args.fg_roc),
    #         Text("C", 25, 20, size=12, weight='bold')
    #         ).move(0, 120),
    #        Panel(
    #         Image(120, 120, args.clinvar_roc),
    #         Text("D", 25, 20, size=12, weight='bold')
    #         ).move(120, 120),
    #        ).save(args.out)

if __name__ == "__main__":
    desc = 'composse grant figure'
    parser = argparse.ArgumentParser(description=desc)
    argLs = ('gene_counts', 'mpc_dist', 'fg_roc', 'clinvar_roc', 'out')
    for param in argLs:
        parser.add_argument(param)
    args = parser.parse_args()
    main(args)

