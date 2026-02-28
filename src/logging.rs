use std::time::{SystemTime, UNIX_EPOCH};

pub fn info(message: impl AsRef<str>) {
    let ts = now_unix_millis();
    let human = format_rfc3339_utc(ts);
    // ANSI green for INFO
    println!("{human}  \x1b[32mINFO\x1b[0m {}", message.as_ref());
}

pub fn warn(message: impl AsRef<str>) {
    let ts = now_unix_millis();
    let human = format_rfc3339_utc(ts);
    // ANSI yellow for WARN
    println!("{human}  \x1b[33mWARN\x1b[0m {}", message.as_ref());
}

fn now_unix_millis() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn format_rfc3339_utc(unix_ms: u128) -> String {
    let secs = (unix_ms / 1000) as i64;
    let ms = (unix_ms % 1000) as u32;

    let days = secs.div_euclid(86_400);
    let sod = secs.rem_euclid(86_400);
    let hour = (sod / 3600) as u32;
    let min = ((sod % 3600) / 60) as u32;
    let sec = (sod % 60) as u32;

    let (year, month, day) = civil_from_days(days);
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{min:02}:{sec:02}.{ms:03}000Z")
}

fn civil_from_days(days_since_unix_epoch: i64) -> (i32, u32, u32) {
    // Howard Hinnant's civil-from-days algorithm; epoch reference is 1970-01-01.
    let z = days_since_unix_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097; // [0, 146096]
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365; // [0, 399]
    let mut y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = mp + if mp < 10 { 3 } else { -9 }; // [1, 12]
    y += if m <= 2 { 1 } else { 0 };
    (y as i32, m as u32, d as u32)
}
