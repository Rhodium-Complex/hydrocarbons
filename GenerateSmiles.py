import numpy as np
import itertools
import re
from openbabel import pybel
from concurrent.futures import ProcessPoolExecutor


class Molecule:
    def __init__(self, bonds):
        self.bonds = bonds
        self.eig = np.around(np.sort(np.linalg.eigvalsh(self.bonds)).real, decimals=1)
        self.hash = tuple(self.eig)
    def __len__(self):
        return len(self.bonds)

def create_C_set(c_num, h_num):
    def generate_combinations(c_num, h_num, current_combination):
        if c_num == 0 and h_num == 0:
            yield np.array(current_combination)
            return
        if c_num < 0 or h_num < 0:
            return

        tmp = 1
        if len(current_combination) > 0:
            tmp = current_combination[0]
            
        # 隣接する水素の数が1、2、3、4のそれぞれのケースを処理
        for i in range(tmp, 5):
            new_combination = [i] + current_combination
            yield from generate_combinations(c_num - 1, h_num - (4 - i), new_combination)

    yield from generate_combinations(c_num, h_num, [])

def create_single_bonds_map(C_type_list): 
    def create_C_map_child(adj_matrix, b, is_unvisitable, start_index = 0):
        unprocessed = np.where(b > 0)[0]  # 未処理の最初のindexを取得。
        if len(unprocessed) == 0:         # 全点が次数を満足している場合、探索終了。
            yield adj_matrix
            return
        i = unprocessed[0]

        # 未処理の点が最後の1つだった場合、探索失敗。
        if i == len(b): return

        for j in range(max(i + 1, start_index), len(b)):
            # 単結合のみの化合物を対象とするため、多重結合はスキップ
            if adj_matrix[j][i] == 1: continue
                
            # 結合先候補に余裕がない場合スキップ
            if b[j] == 0: continue
            
            # 候補が探索済みもしくは未探索のメチン、メチレン、メタンの最初でなければスキップ
            if is_unvisitable[j] and j-i>1: continue

            (tmp_a, tmp_b, tmp_visited) = (adj_matrix.copy(), b.copy(), is_unvisitable.copy())

            (tmp_a[i][j], tmp_a[j][i]) = (1, 1)
            (tmp_b[i], tmp_b[j]) = (tmp_b[i] - 1, tmp_b[j] - 1)
            tmp_visited[j + 1] = False

            yield from create_C_map_child(
                tmp_a, 
                tmp_b, 
                tmp_visited, 
                j if tmp_b[i] else 0)
    
    C_len = len(C_type_list)
    visited = (
        [False]
        + [C_type_list[i + 1] == C_type_list[i] for i in range(C_len - 1)]
        + [True]
    )
    visited[1] = False
    
    yield from create_C_map_child(
        np.full((C_len, C_len), 0, dtype="i4"),
        C_type_list,
        visited
    )

def is_connected_graph(adj_matrix):
    """
    Check if a graph is connected.
    """
    visited = np.zeros(len(adj_matrix), dtype=bool)
    visited[0] = True
    for _ in range(len(adj_matrix) - 1):
        visited |= np.dot(adj_matrix, visited).astype(bool)
        if visited.all():
            return True
    return False

def Morgan(m_mat):
    connect_map = m_mat#np.array(m_mat > 0)
    flag = 0
    rst = 4 ** np.sum(connect_map, axis = 1)
    tmp = np.unique(rst)
    while len(tmp) > flag:
        flag = len(tmp)
        for i in range(len(tmp)):
            rst[rst == tmp[i]] = 4**(flag-i)
        rst = connect_map @ rst
        tmp = np.unique(rst)
    return rst

def canonicalize(m_mat: np.ndarray) -> np.ndarray:
    is_connected = (m_mat > 0)
    cano_num = Morgan(m_mat)

    def gen():
        for i in range(11)[::-1]:
            yield i 
    gen = gen()

    tmp = np.zeros(len(cano_num))
    m_mat2 = (m_mat > 0) * cano_num + m_mat
    tmp[cano_num.argmax()] = 11

    stack = [i for i in m_mat2[cano_num.argmax()].argsort() if m_mat2[cano_num.argmax()][i]]
    for i in stack:
        if tmp[i] == 0 :
            tmp[i] = next(gen)
            stack += [j for j in m_mat2[i].argsort() if m_mat2[i][j]]
    cano_num = tmp * 10

    a   = cano_num.argmax()
    rst = [a]
    cano_num[a] = 0
    asd = is_connected[a] * cano_num + m_mat[a]

    #ラベリングが一番大きい枝を構成
    while max(asd) > 10:
        a = asd.argmax()
        rst = [a] + rst
        cano_num[a] = 0
        asd = is_connected[a] * cano_num + m_mat[a]

    #ラベリング2番目、及び小枝を構成
    counter = len(rst) - 1
    asd = is_connected[rst[counter]] * cano_num + m_mat[rst[counter]]
    if max(asd)<=10:
        counter = 0
    while max(cano_num) > 10:
        while max(asd) <= 10:
            if counter >= len(rst):
                break
            asd = is_connected[rst[counter]] * cano_num + m_mat[rst[counter]]
            counter += 1
        a = asd.argmax()
        rst = rst + [a]
        cano_num[a] = 0
        asd = is_connected[a] * cano_num
        counter += 1
        
        if max(asd) == 0:
            counter = 0
    return np.array(m_mat[:, rst][rst,:])

def unique_mols(m_mat: list):
    rst = {}
    for x in m_mat:
        fing_m = x.hash
        fing = Morgan(x.bonds)
        if len(np.unique(fing))==1: isHZreo=True
        fing = [list(np.where(fing == i)[0]) for i in np.unique(fing)[::-1]]
        
        tmpNull = sum(next(itertools.product(*[ itertools.permutations(j) for j in fing])),())
        if fing_m not in rst:
            rst[fing_m] = [x.bonds[:,tmpNull][tmpNull,:]]
            yield x
            continue
        
        for i in itertools.product(*[ itertools.permutations(j) for j in fing]):
            i = sum(i,())
            tmp = x.bonds[:,i][i,:]
            if any(np.array_equal(tmp, j) for j in rst[fing_m]): break
        else:
            rst[fing_m].append(x.bonds[:,tmpNull][tmpNull,:])
            yield x
    del rst

def mat2smiles(mat):
    mat = mat.bonds
    char    = ["","","=","#"]
    rst     = [""] * len(mat[0])
    counter = 1

    for i in range(len(mat[0]))[::-1]:
        flag = False
        rst[i]="C" + rst[i]
        for j in range(i)[::-1]:
            bond = mat[i][j]
            if bond == 0:
                continue
            if flag:
                rst[i] += char[bond] + "{"+str(counter)+"}"
                rst[j] = "{"+str(counter)+"}"+rst[j]
                counter += 1
            else:
                rst[i] = char[bond] + rst[i]
            flag = True
        tmp = ""
        flag = False
        for j in range(i, len(mat[0]))[::-1]:
            if mat[i][j] == 0:
                continue
            if rst[j] == "":
                continue
            tmp = "(" + rst[j] + ")" + tmp if flag else rst[j]
            rst[j] = ""
            flag = True
        rst[i] += tmp
    rst = rst[0]
    
    counter = 1
    def sanit(num):
        if num < 10:return str(num)
        return "%"+str(num)
    
    while rst.find("{") > 0:
        i = rst.find("{")
        l = rst.find("}") 
        
        pattern =    '((?:[=#]{0,1}\{\d})*)(={0,1})\\'+rst[i:l]+'}((?:[=#]*\d)*)'
        replacement = r'\g<3>\g<2>' + sanit(counter) + r'\g<1>'
        rst = re.sub(pattern, replacement, rst)
        counter = counter + 1
    return rst


def unique_dehydro_mols(dehydro_stream)-> list:
    def dehydrogenase(mol):
        MAX_BOND_ORDER = 3
        
        for i in range(len(mol)):
            if sum(mol.bonds[i]) > MAX_BOND_ORDER: continue
            for j in range(i + 1, len(mol)):
                if mol.bonds[i][j] == 0: continue
                if mol.bonds[i][j] == MAX_BOND_ORDER: continue
                if sum(mol.bonds[j]) > MAX_BOND_ORDER: continue
                
                tmp = mol.bonds.copy()
                tmp[i][j] +=1
                tmp[j][i] +=1
                yield Molecule(tmp)
                
    def generate_molecules(j):
        for row in j:
            yield from dehydrogenase(row)
    
    return list(unique_mols(generate_molecules(dehydro_stream)))

def build_structure(j):
    candidate_mols = [canonicalize(i) for i in create_single_bonds_map(j) if is_connected_graph(i)]
    unique_mols_list = unique_mols(Molecule(i) for i in np.unique(candidate_mols, axis=0))
    return [[Moleclue] for Moleclue in unique_mols_list]

#C8までで2分半
#C9までで17分ぐらい
#def test():

tmp2 = []
sum_rst = ['N#N', #C1H0
        'N#N', #C1H2
        'N#N', #C1H4
        'C',]
for C in range(2,11):
    tmp = []
    print(C, ':', end = '')
    for H in range(0, C * 2 + 3, 2)[::-1]:
        sum_rst.append('N#N')
        if __name__ == '__main__':
            with ProcessPoolExecutor() as executor:
                future = executor.map(unique_dehydro_mols, tmp)
            tmp = list(future)
            print('.', end = '')
            with ProcessPoolExecutor() as executor: 
                future = executor.map(build_structure, create_C_set(C,H))
            print('.', end = '')
            for i in future:
                tmp += i

        tmp2 = (x for row in tmp for x in row)
        future = map(mat2smiles, tmp2)
        sum_rst_before = len(sum_rst)
        sum_rst += list(future)
        print(len(sum_rst)-sum_rst_before, end = ' ')

    print('')
