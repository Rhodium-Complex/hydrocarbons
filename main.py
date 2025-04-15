from openbabel import pybel
from concurrent.futures import ProcessPoolExecutor
import structure_generator, molecule_transformations, converter
import itertools


#C8までで2分半
#C9までで17分ぐらい
#def test(): 

def main ():
    with ProcessPoolExecutor() as dehydro_executor, ProcessPoolExecutor() as structure_executor:
        # initialization: 最初の4つは同じ
        sum_rst = ['N#N', #C1H0
                'N#N', #C1H2
                'N#N', #C1H4
                'C',]
        for C in range(2,8):
            molecule_structures = []
            print(C, ':', end = '')
            for H in range(0, C * 2 + 3, 2)[::-1]:
                # ラベルを追加
                sum_rst.append('N#N')
                
                future = dehydro_executor.map(
                    molecule_transformations.unique_dehydro_mols, 
                    molecule_structures
                )
                molecule_structures = list(future)
                print('.', end = '')
                    
                future = structure_executor.map(
                    structure_generator.build_structure, 
                    structure_generator.build_carbon_hydrogen_combination(C,H)
                )
                print('.', end = '')
                for i in future:
                    molecule_structures += i

                flattened_structures = itertools.chain.from_iterable(molecule_structures)
                future = map(converter.mat2smiles, flattened_structures)
                sum_rst_before = len(sum_rst)
                sum_rst += list(future)
                print(len(sum_rst)-sum_rst_before, end = ' ')

            print('')
    return sum_rst

# ファイル全体のメイン部分をここで定義
if __name__ == '__main__':
    sum_rst = main()


# 2 :1 1 1 0 
# 3 :1 2 3 2 1 
# 4 :2 5 9 11 7 3 
# 5 :3 10 26 40 40 21 6 
# 6 :5 25 77 159 217 185 85 19 
# 7 :9 56 222 574 1029 1229 920 356 50 
# 8 :18 139 652 2069 4656 7396 7950 5289 1804 204

#8/28 21:50-9/1 22:33
#2 :...1 ...1 ...1 ...
#3 :...1 ...2 ...3 ...2 ...1 
#4 :...2 ...5 ...9 ...11 ...7 ...3 
#5 :...3 ...10 ...26 ...40 ...40 ...21 ...6 
#6 :...5 ...25 ...77 ...159 ...217 ...185 ...85 ...19 
#7 :...9 ...56 ...222 ...575 ...1031 ...1230 ...920 ...356 ...50 
#8 :...18 ...139 ...654 ...2082 ...4679 ...7437 ...7982 ...5308 ...1804 ...204 
#9 :...35 ...338 ...1902 ...7244 ...19983 ...40139 ...57771 ...56437 ...33860 ...10064 ...832 
#10 :...75 ...852 ...5568 ...24938 ...81909 ...201578 ..369067 ..488125 ..439373 ..241297 ..64352 .

# 作成中のため強制終了
exit()

import pickle
with open('mols', 'wb') as f:
  pickle.dump(sum_rst , f)



print(sum_rst[27073])
sum_rst[27073] = "CC"
print(sum_rst[27297])
sum_rst[27297] = "CC"
print(sum_rst[31837])
sum_rst[31837] = "CC"
print(sum_rst[31995])
sum_rst[31995] = "CC"

svg_group_size = 190

for i,mols in enumerate([sum_rst[idx:idx + svg_group_size] for idx in range(0,len(sum_rst), svg_group_size)]):
    tmp = ['<g transform="translate({0},{1})">\n'.format(ji%10*100, int(ji/10)*100)+ pybel.readstring("smi", j
        ).write('svg', opt={'x':True}).split(
        '<g transform="translate(0,0)">\n')[1].split('</g>')[0].replace('opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"', '')+'</g>' for (ji,j) in enumerate(mols)]
    with open('./image/{0}.svg'.format(i), mode='w') as f:
        f.write("""<?xml version="1.0"?>
<svg version="1.1" id="topsvg"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
xmlns:cml="http://www.xml-cml.org/schema" x="0" y="0" width="1000px" height="1900px" viewBox="0 0 1000 1900">
<title> - Open Babel Depiction</title>
<rect x="0" y="0" width="100" height="100" fill="white"/>"""+"\n".join(tmp)+'</svg>')

for i in range(26790,100000):
    print(i,pybel.readstring("smi", sum_rst[i]))
    pybel.readstring("smi", sum_rst[i]).write('svg')

