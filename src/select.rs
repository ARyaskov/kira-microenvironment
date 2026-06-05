//! Order-statistic helpers (median, quantile, trimmed mean) on `&mut [f32]`.
//!
//! Implemented with `select_nth_unstable_by` (O(n) expected) instead of full sort,
//! which is the dominant cost in Stage1/Stage5 aggregation. All routines remain
//! deterministic because `f32::total_cmp` provides a total order including NaN.

#[inline]
fn cmp_f32(a: &f32, b: &f32) -> std::cmp::Ordering {
    a.total_cmp(b)
}

/// Median of `values`. Mutates in place (reorders) but does not fully sort.
pub fn median_in_place(values: &mut [f32]) -> f32 {
    let n = values.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return values[0];
    }
    if n % 2 == 1 {
        let mid = n / 2;
        let (_, m, _) = values.select_nth_unstable_by(mid, cmp_f32);
        *m
    } else {
        let hi = n / 2;
        // Upper of the two middle elements.
        let (lo_half, m_hi, _) = values.select_nth_unstable_by(hi, cmp_f32);
        let upper = *m_hi;
        // Max of the lower half is the lower median element.
        // (After `select_nth_unstable`, lo_half holds the n/2 smallest in unspecified order.)
        let lower = lo_half
            .iter()
            .copied()
            .reduce(|a, b| if a.total_cmp(&b).is_ge() { a } else { b })
            .unwrap_or(upper);
        (lower + upper) * 0.5
    }
}

/// Trimmed mean. `trim` is clamped to `[0, 0.499999]`.
///
/// Partitions twice via `select_nth_unstable` (O(n)) instead of full sort.
pub fn trimmed_mean_in_place(values: &mut [f32], trim: f32) -> f32 {
    let n = values.len();
    if n == 0 {
        return 0.0;
    }
    let t = trim.clamp(0.0, 0.499_999);
    let k = ((n as f32) * t).floor() as usize;
    if k * 2 >= n {
        return values.iter().copied().sum::<f32>() / (n as f32);
    }

    // Partition so that elements with rank < k are in the left segment.
    let (_, _, right_inclusive) = values.select_nth_unstable_by(k, cmp_f32);
    // `right_inclusive` covers indices [k+1, n), but we want the slice [k, n-k).
    // Operate on the original `values` after the partition:
    let n_keep_right = n - k;
    if n_keep_right <= k + 1 {
        // Defensive; shouldn't happen because k*2 < n.
        return values[k];
    }
    // Re-partition the right segment (which now starts at index k) so that the
    // last `k` elements are the largest. The remaining middle `n - 2k` is what we average.
    let mid_len = n - 2 * k;
    if mid_len == 0 {
        return values[k];
    }
    // Within values[k..], rank mid_len-1 separates the trimmed-mean window from the upper trim.
    let _ = right_inclusive; // borrow released
    let right = &mut values[k..];
    let _ = right.select_nth_unstable_by(mid_len - 1, cmp_f32);
    let window = &values[k..k + mid_len];
    let sum: f32 = window.iter().copied().sum();
    sum / (mid_len as f32)
}

/// `p`-quantile (linear interpolation NOT applied — uses ceil-rank to match prior behavior).
///
/// Rank: `rank = ceil(p * n) - 1`, clamped to `[0, n-1]`.
pub fn quantile_in_place(values: &mut [f32], p: f32) -> f32 {
    let n = values.len();
    if n == 0 {
        return 0.0;
    }
    let pp = p.clamp(0.0, 1.0);
    let rank = ((pp * n as f32).ceil() as usize)
        .saturating_sub(1)
        .min(n - 1);
    let (_, v, _) = values.select_nth_unstable_by(rank, cmp_f32);
    *v
}
