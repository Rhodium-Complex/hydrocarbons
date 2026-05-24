"""Module for converting molecule bond matrices into SMILES-like strings."""
def mat2smiles(mat):
    """Convert a molecule bond matrix into a compact SMILES-like string."""
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

    while "{" in smiles_parts:
        i = smiles_parts.find("{")
        l = smiles_parts.find("}")
        ring_token = smiles_parts[i:l + 1]

        # 正規表現を使用してラベルを置換
        # 例: {1} を適切な SMILES 記法に変換
        smiles_parts = smiles_parts.replace(ring_token, sanit(counter), 2)
        counter = counter + 1
    return smiles_parts
