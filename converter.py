import re

def mat2smiles(mat):
    bond_matrix = mat.bonds
    bond_symbols = ["", "", "=", "#"]
    smiles_parts = [""] * len(bond_matrix[0])
    counter = 1

    for i in range(len(bond_matrix[0]))[::-1]:
        flag = False
        smiles_parts[i] = "C" + smiles_parts[i]
        for j in range(i)[::-1]:
            bond = bond_matrix[i][j]
            if bond == 0:
                continue
            if flag:
                smiles_parts[i] += bond_symbols[bond] + "{" + str(counter) + "}"
                smiles_parts[j] = "{" + str(counter) + "}" + smiles_parts[j]
                counter += 1
            else:
                smiles_parts[i] = bond_symbols[bond] + smiles_parts[i]
            flag = True
        tmp = ""
        flag = False
        for j in range(i, len(bond_matrix[0]))[::-1]:
            if bond_matrix[i][j] == 0:
                continue
            if smiles_parts[j] == "":
                continue
            tmp = "(" + smiles_parts[j] + ")" + tmp if flag else smiles_parts[j]
            smiles_parts[j] = ""
            flag = True
        smiles_parts[i] += tmp
    smiles_parts = smiles_parts[0]
    
    counter = 1
    def sanit(num):
        if num < 10:
            return str(num)
        return "%" + str(num)
    
    while smiles_parts.find("{") > 0:
        i = smiles_parts.find("{")
        l = smiles_parts.find("}") 
        
        # 正規表現を使用してラベルを置換
        # 例: {1} を適切な SMILES 記法に変換
        pattern = r'((?:[=#]{0,1}\{\d})*)(={0,1})\\' + smiles_parts[i:l] + r'}((?:[=#]*\d)*)'
        replacement = r'\g<3>\g<2>' + sanit(counter) + r'\g<1>'
        smiles_parts = re.sub(pattern, replacement, smiles_parts)
        counter = counter + 1
    return smiles_parts
