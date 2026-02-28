pub fn median_in_place(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) * 0.5
    }
}

pub fn trimmed_mean_in_place(values: &mut [f32], trim: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let t = trim.clamp(0.0, 0.499_999);
    let n = values.len();
    let k = ((n as f32) * t).floor() as usize;
    if k * 2 >= n {
        return values.iter().copied().sum::<f32>() / (n as f32);
    }
    let slice = &values[k..(n - k)];
    slice.iter().copied().sum::<f32>() / (slice.len() as f32)
}

pub fn quantile_in_place(values: &mut [f32], p: f32) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let pp = p.clamp(0.0, 1.0);
    let n = values.len();
    let rank = ((pp * n as f32).ceil() as usize)
        .saturating_sub(1)
        .min(n - 1);
    values[rank]
}
