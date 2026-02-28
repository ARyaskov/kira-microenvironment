#[derive(Debug, Clone, Copy, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum SimdLevel {
    Avx2,
    Neon,
    Scalar,
}

impl SimdLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Avx2 => "avx2",
            Self::Neon => "neon",
            Self::Scalar => "scalar",
        }
    }
}

pub fn detect_simd_level() -> SimdLevel {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            return SimdLevel::Avx2;
        }
        return SimdLevel::Scalar;
    }
    #[cfg(target_arch = "aarch64")]
    {
        return SimdLevel::Neon;
    }
    #[allow(unreachable_code)]
    SimdLevel::Scalar
}
