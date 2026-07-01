# 構造マップ

この文書は、炭化水素生成プログラムの大まかな呼び出し関係をまとめたものです。
処理の流れが分からなくなったときは、まずこのファイルを見てください。

生成 pipeline を変更したときは、この文書も同じタイミングで更新してください。

## 入口

コマンド実行時の入口は `main.py` です。

```text
python main.py
  -> parse_args()
  -> main(...)
```

`main.main()` は、指定されたオプションに応じて3つの経路に分かれます。

```text
main(...)
  --export pdf/svg の場合:
    -> generation_pipeline.run_export_smiles_groups(...)
    -> export_structures_pdf.export_formula_smiles_pdf(...)
       または export_structures_svg.export_formula_smiles_svg_pages(...)

  include_smiles=True の場合:
    -> generation_pipeline.run_generation_smiles_groups(...)

  include_smiles=False の場合:
    -> generation_pipeline.run_generation(...)
```

## 生成処理の中心

`run_generation()` と `run_generation_smiles_groups()` は公開 API としては分かれていますが、
内部ではどちらも `_run_generation_pipeline()` を使います。

```text
generation_pipeline.run_generation(...)
  -> _run_generation_pipeline(include_smiles=False, ...)

generation_pipeline.run_generation_smiles_groups(...)
  -> _run_generation_pipeline(include_smiles=True, ...)
```

`_run_generation_pipeline()` は共通の生成 loop を持ち、SMILES を作るかどうかだけを切り替えます。

```text
generation_pipeline._run_generation_pipeline(...)
  -> _iter_formula_structure_steps(...)
     -> molecule_transformations.unique_dehydro_mols(...)
     -> structure_generator.build_carbon_hydrogen_combination(...)
     -> structure_generator.build_structure(...)
```

`include_smiles=True` の場合は、構造生成の後に SMILES 変換が追加されます。

```text
generation_pipeline._run_generation_pipeline(include_smiles=True, ...)
  -> _structure_smiles(...)
     -> converter.mat2smiles_variants(...)
        -> converter._render_smiles(...)
        -> stereochemistry.analyze_ez(...)      [include_stereo=True のときだけ]
        -> converter.mat2stereo_smiles(...)     [include_stereo=True のときだけ]
```

`_iter_formula_structure_steps()` は、炭素数ごとに、水素数を多い方から少ない方へ順に処理します。
各水素数ステップでは、まず前段で得た構造を脱水素化し、その後で現在の分子式に対応する単結合骨格を追加します。

```text
for carbon_count in min_carbon..max_carbon:
  current_carbon_structures = []
  for hydrogen_count in 多い方から少ない方へ:
    current_carbon_structures = 前の構造を脱水素化する
    current_carbon_structures += 現在の C/H 分子式の単結合骨格を作る
    1つの分子式ステップとして返す
```

## 単結合骨格の生成

まず、分子式から「各炭素が他の炭素と作る単結合の本数」の候補を作ります。

```text
structure_generator.build_carbon_hydrogen_combination(c_num, h_num)
  -> single_bond_degrees 配列を返す
```

`single_bond_degrees` の各値は、初期の単結合骨格で、その炭素が必要とする C-C 単結合数です。

次に、その次数列を満たす単結合隣接行列を作ります。

```text
structure_generator.build_structure(single_bond_degrees)
  -> create_single_bonds_map(single_bond_degrees)
     -> recursive generate_child_nodes(...)
        -> graph_utils.is_connected_graph(...)
  -> graph_utils.canonicalize(...)
  -> np.unique(...)
  -> deduplication.unique_mols(...)
     -> graph_utils.morgan(...)
     -> isomorphism.has_permutation_match(...)
  -> [[Molecule], [Molecule], ...] を返す
```

`create_single_bonds_map()` が作るのは単結合だけの骨格です。
二重結合・三重結合は、この後の脱水素化で作られます。

## 脱水素化

脱水素化は `molecule_transformations.unique_dehydro_mols()` が担当します。

```text
molecule_transformations.unique_dehydro_mols(structures)
  -> generate_candidates(...)
     -> 元の結合行列をコピーする
     -> 既存の C-C 結合を1つ選んで結合次数を +1 する
     -> 結合がない場所、三重結合、価数超過は除外する
  -> deduplication.unique_mols(...)
```

この処理により、`C_n H_m` の構造から `C_n H_{m-2}` の候補が作られます。
二重結合と三重結合はここで発生します。

## SMILES 変換

通常の SMILES 出力は次の流れです。

```text
generation_pipeline._structure_smiles(...)
  -> converter.mat2smiles_variants(...)
     -> converter._render_smiles(...)
```

E/Z 立体を含める場合は、追加で stereochemistry を使います。

```text
converter.mat2smiles_variants(..., include_stereo=True)
  -> stereochemistry.analyze_ez(...)
  -> converter.mat2stereo_smiles(...)
```

ここでの E/Z は正式な CIP 命名ではなく、二重結合の二つの相対配置を
列挙するための内部ラベルです。環状構造では、立体制約を先に解いてから
全域木と環閉鎖辺の両方へ `/`・`\\` を付与します。制約が矛盾する、標準
SMILES で表現できない、または同じ立体文字列へ重複した候補は、列挙位置を
維持するため立体情報なしの SMILES にフォールバックします。

## PDF / SVG 出力

PDF 出力:

```text
main(..., export_format="pdf")
  -> generation_pipeline.run_export_smiles_groups(...)
  -> export_structures_pdf.export_formula_smiles_pdf(...)
     -> structure_export_layout.build_structure_cells(...)
     -> RDKit で SMILES を分子図に変換
     -> ReportLab で PDF を出力
```

SVG 出力:

```text
main(..., export_format="svg")
  -> generation_pipeline.run_export_smiles_groups(...)
  -> export_structures_svg.export_formula_smiles_svg_pages(...)
     -> structure_export_layout.build_structure_cells(...)
     -> RDKit で SMILES を分子図に変換
     -> SVG ページを出力
```

## モジュールの役割

- `main.py`: CLI 引数の解釈と実行モードの選択。
- `generation_pipeline.py`: 全体の orchestration、計時、並列処理、出力グループ化。
- `structure_generator.py`: 分子式から次数列を作り、連結な単結合骨格を生成する。
- `molecule_transformations.py`: 既存結合の結合次数を上げて脱水素化する。
- `deduplication.py`: 構造重複を除去する。
- `graph_utils.py`: 連結判定、Morgan ラベル、canonicalize。
- `isomorphism.py`: 同型な結合行列の permutation 判定。
- `molecule.py`: 結合行列を持つ `Molecule` と fingerprint。
- `converter.py`: 結合行列から SMILES 風文字列を作る。
- `stereochemistry.py`: 形式的な E/Z assignment を展開する。
- `export_structures_pdf.py`: PDF 出力。
- `export_structures_svg.py`: SVG 出力。
- `structure_export_layout.py`: B5 グリッド配置と RDKit warning 周りの共通処理。

## 更新ルール

次のどれかを変更した場合は、この文書も更新してください。

- `main.py` の公開入口や CLI モード。
- 脱水素化、単結合骨格生成、重複排除、SMILES 変換の順序。
- `single_bond_degrees` の意味。
- `create_single_bonds_map()` が返すものの条件。
- deduplication、canonicalization、stereochemistry、export layout の担当モジュール。

処理の流れに関する説明は、各モジュールに散らすより、このファイルに短い call map として追加してください。
