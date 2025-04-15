from openbabel import pybel
from concurrent.futures import ProcessPoolExecutor
import structure_generator, molecule_transformations, converter
import itertools
import os # For creating output directory

# --- Configuration Constants ---
MIN_CARBON = 2
MAX_CARBON = 7 
SVG_GROUP_SIZE = 190
SVG_OUTPUT_DIR = './image/'
# ---

#C8までで2分半
#C9までで17分ぐらい
#def test(): 

def main ():
    # dehydro_executor: 脱水素化反応による分子構造生成を担当
    # structure_executor: 炭素と水素の組み合わせから基本的な構造を生成を担当
    with ProcessPoolExecutor() as dehydro_executor, ProcessPoolExecutor() as structure_executor:
        # initialization: 全ての分子構造のSMILES文字列を格納するリスト
        # 初期値は C=1 のケースに対応。'N#N' は該当する安定構造がない場合のプレースホルダー。
        # C1H0, C1H2, C1H4 は安定構造なし、C1 (Hなし) は 'C' (メタンラジカル等だがここでは単純に炭素原子として)
        all_smiles_results = ['N#N', # C1H0 (該当なし)
                             'N#N', # C1H2 (該当なし)
                             'N#N', # C1H4 (該当なし)
                             'C',   # C1 (メタン CH4 だが、ここでは炭素骨格のみ考慮か？要確認)
                             ]
        # 炭素数 MIN_CARBON から MAX_CARBON までループ
        for C in range(MIN_CARBON, MAX_CARBON + 1):
            # 現在の炭素数 C における分子構造候補 (SMILES生成前の中間表現) を一時的に格納するリスト
            current_carbon_structures = []
            print(C, ':', end = '')
            for H in range(0, C * 2 + 3, 2)[::-1]:
                # Hに対応するインデックスにプレースホルダーを追加 (後で上書きされる)
                all_smiles_results.append('N#N')

                # 1. 脱水素化による構造生成 (既存の構造リストから新しい構造を生成)
                future_dehydro = dehydro_executor.map(
                    molecule_transformations.unique_dehydro_mols,
                    current_carbon_structures # 前のHステップ、または初期の空リストを使用
                )
                # mapの結果をリストに変換して更新
                current_carbon_structures = list(future_dehydro)
                print('.', end = '')

                # 2. C/H組み合わせからの基本構造生成
                future_structure = structure_executor.map(
                    structure_generator.build_structure,
                    structure_generator.build_carbon_hydrogen_combination(C,H)
                )
                print('.', end = '')
                # structure_executorの結果を既存のリストに追加
                for structures in future_structure:
                    current_carbon_structures += structures

                # 3. SMILESへの変換と結果リストへの追加
                #    current_carbon_structures はリストのリストになっている可能性があるため、フラット化
                flattened_structures = itertools.chain.from_iterable(current_carbon_structures)
                future_smiles = map(converter.mat2smiles, flattened_structures)
                # この C, H の組み合わせで生成されたSMILES数をカウントするために、追加前のリスト長を記録
                results_before_adding = len(all_smiles_results)
                # 生成されたSMILES文字列をリストに追加
                all_smiles_results += list(future)
                # 追加されたSMILESの数を表示
                print(len(all_smiles_results) - results_before_adding, end=' ')

            print('')
    return all_smiles_results

# ファイル全体のメイン部分をここで定義
if __name__ == '__main__':
    all_smiles_results = main()


# --- Example Output Counts (C: H=...) ---
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
exit() # Still keeping the exit() as per user denial

import pickle
# 生成されたSMILESリストをpickleファイルに保存
with open('mols.pkl', 'wb') as f: # Changed extension to .pkl for clarity
  pickle.dump(all_smiles_results , f)


# # --- Debugging/Testing Section (Commented out) ---
# # This section appears to be for debugging specific indices or testing modifications.
# # It's commented out for regular execution.
# print(all_smiles_results[27073])
# all_smiles_results[27073] = "CC" # Example modification
# print(all_smiles_results[27297])
# all_smiles_results[27297] = "CC" # Example modification
# print(all_smiles_results[31837])
# all_smiles_results[31837] = "CC" # Example modification
# print(all_smiles_results[31995])
# all_smiles_results[31995] = "CC" # Example modification
# # --- End Debugging/Testing Section ---


# --- SVG Generation ---
# Ensure the output directory exists
os.makedirs(SVG_OUTPUT_DIR, exist_ok=True)

# SVG生成ループ: all_smiles_results を SVG_GROUP_SIZE ごとに分割して処理
for i,mols in enumerate([all_smiles_results[idx:idx + SVG_GROUP_SIZE] for idx in range(0,len(all_smiles_results), SVG_GROUP_SIZE)]):
    svg_parts = []
    for ji, j_smiles in enumerate(mols):
        try:
            # SMILES文字列からOpen Babel分子オブジェクトを生成し、SVGデータを取得
            mol = pybel.readstring("smi", j_smiles)
            # TODO: Clarify the meaning of opt={'x':True} and the string replacement below
            svg_content = mol.write('svg', opt={'x':True}).split( # opt={'x':True} の意味は要確認
                '<g transform="translate(0,0)">\n')[1].split('</g>')[0].replace('opacity="1.0" stroke="rgb(0,0,0)"  stroke-width="2.0"', '') # この置換の理由も要確認
            # グループ化してリストに追加
            svg_parts.append('<g transform="translate({0},{1})">\n'.format(ji%10*100, int(ji/10)*100) + svg_content + '</g>')
        except Exception as e:
            print(f"Warning: Failed to generate SVG for SMILES '{j_smiles}' at index {i * SVG_GROUP_SIZE + ji}. Error: {e}")
            # エラーが発生した場合、空のグループを追加するか、何も追加しないかを選択できます。
            # ここでは何も追加しないことにします。
            pass
    tmp = svg_parts # Assign the collected parts to tmp
    # SVGファイルへの書き込み
    output_svg_path = os.path.join(SVG_OUTPUT_DIR, f'{i}.svg')
    with open(output_svg_path, mode='w') as f:
        f.write("""<?xml version="1.0"?>
<svg version="1.1" id="topsvg"
xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
xmlns:cml="http://www.xml-cml.org/schema" x="0" y="0" width="1000px" height="1900px" viewBox="0 0 1000 1900">
<title> - Open Babel Depiction</title>
<rect x="0" y="0" width="100" height="100" fill="white"/>"""+"\n".join(tmp)+'</svg>')

