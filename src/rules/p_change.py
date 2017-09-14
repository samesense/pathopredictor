from Bio.Data import IUPACData

def convert_protein_change(protein_change):
    # can't handle termination yet
    if 'p.' in protein_change and not 'Ter' in protein_change:
        # p.Ile199Val
        # p.Arg151*
        # p.Val276fs
        p1 = protein_change.split('.')[1][:3]
        print(protein_change)
        if '?' in protein_change:
            c2 = ''
            protein_pos = protein_change.split('.')[1][3:-1]
            c1 = IUPACData.protein_letters_3to1[p1]
        elif 'del' in protein_change:
            if '_' in protein_change:
                c2 = ''
                protein_pos = protein_change.split('.')[1].split('_')[0][3:]
                c1 = IUPACData.protein_letters_3to1[p1]
            else:
                c2 = ''
                protein_pos = protein_change.split('.')[1][3:-3]
                c1 = IUPACData.protein_letters_3to1[p1]
        elif 'fs' in protein_change:
            c2 = ''
            protein_pos = protein_change.split('.')[1][3:-2]
            c1 = IUPACData.protein_letters_3to1[p1]
        elif '*' in protein_change:
            c2 = ''
            protein_pos = protein_change.split('.')[1][3:-1]
            c1 = IUPACData.protein_letters_3to1[p1]
        else:    
            p2 = protein_change.split('.')[1][-3:]
            protein_pos = protein_change.split('.')[1][3:-3]
            c1 = IUPACData.protein_letters_3to1[p1]
            c2 = IUPACData.protein_letters_3to1[p2]
        return c1 + protein_pos + c2
    else:
        return 'NA'
