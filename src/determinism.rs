use std::cmp::Ordering;

pub fn stable_sort_by<T, F>(items: &mut [T], cmp: F)
where
    F: FnMut(&T, &T) -> Ordering,
{
    items.sort_by(cmp);
}

pub fn cmp_opt_str(a: &Option<String>, b: &Option<String>) -> Ordering {
    match (a, b) {
        (Some(x), Some(y)) => x.cmp(y),
        (None, Some(_)) => Ordering::Less,
        (Some(_), None) => Ordering::Greater,
        (None, None) => Ordering::Equal,
    }
}
