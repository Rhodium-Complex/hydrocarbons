use pyo3::prelude::*;

fn permute_group(group: &[usize], used: &mut [bool], current: &mut Vec<usize>, output: &mut Vec<Vec<usize>>) {
    if current.len() == group.len() {
        output.push(current.clone());
        return;
    }

    for index in 0..group.len() {
        if used[index] {
            continue;
        }
        used[index] = true;
        current.push(group[index]);
        permute_group(group, used, current, output);
        current.pop();
        used[index] = false;
    }
}

fn group_permutations(group: &[usize]) -> Vec<Vec<usize>> {
    let mut output = Vec::new();
    let mut used = vec![false; group.len()];
    let mut current = Vec::with_capacity(group.len());
    permute_group(group, &mut used, &mut current, &mut output);
    output
}

fn matching_prefix_record_indices(
    matrix: &[u8],
    n: usize,
    permutation: &[usize],
    representatives: &[Vec<u8>],
    candidate_indices: &[usize],
) -> Vec<usize> {
    let mut matches = Vec::new();
    for &representative_index in candidate_indices {
        let representative = &representatives[representative_index];
        let mut is_match = true;
        'rows: for (row_position, &row_index) in permutation.iter().enumerate() {
            let matrix_row_offset = row_index * n;
            let representative_row_offset = row_position * n;
            for (column_position, &column_index) in permutation.iter().enumerate() {
                if matrix[matrix_row_offset + column_index]
                    != representative[representative_row_offset + column_position]
                {
                    is_match = false;
                    break 'rows;
                }
            }
        }
        if is_match {
            matches.push(representative_index);
        }
    }
    matches
}

fn search_product(
    matrix: &[u8],
    n: usize,
    groups: &[Vec<Vec<usize>>],
    group_index: usize,
    permutation: &mut Vec<usize>,
    representatives: &[Vec<u8>],
    candidate_indices: &[usize],
) -> bool {
    if group_index == groups.len() {
        return !candidate_indices.is_empty();
    }

    for group_permutation in &groups[group_index] {
        let original_len = permutation.len();
        permutation.extend(group_permutation);
        let next_candidate_indices =
            matching_prefix_record_indices(matrix, n, permutation, representatives, candidate_indices);
        if !next_candidate_indices.is_empty()
            && search_product(
                matrix,
                n,
                groups,
                group_index + 1,
                permutation,
                representatives,
                &next_candidate_indices,
            )
        {
            return true;
        }
        permutation.truncate(original_len);
    }
    false
}

#[pyfunction]
fn has_permutation_match(
    matrix: Vec<u8>,
    n: usize,
    bond_fingerprints: Vec<Vec<usize>>,
    representative_keys: Vec<Vec<u8>>,
) -> PyResult<bool> {
    if matrix.len() != n * n {
        return Err(pyo3::exceptions::PyValueError::new_err("matrix length does not match n*n"));
    }

    let representatives = representative_keys;
    if representatives.is_empty() {
        return Ok(false);
    }

    let groups: Vec<Vec<Vec<usize>>> = bond_fingerprints
        .iter()
        .map(|group| group_permutations(group))
        .collect();
    let mut permutation = Vec::with_capacity(n);
    let candidate_indices = (0..representatives.len()).collect::<Vec<_>>();
    Ok(search_product(
        &matrix,
        n,
        &groups,
        0,
        &mut permutation,
        &representatives,
        &candidate_indices,
    ))
}

#[pymodule]
fn hydrocarbon_rust(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(has_permutation_match, module)?)?;
    Ok(())
}
