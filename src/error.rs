use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    InputMissing,
    TsvHeader,
    TsvParse,
    Path,
    InvalidArgument,
    ResourcesIncomplete,
    GroupsFormat,
    SecretionInputMissing,
    RegimeMapMissing,
    RegimeMapParse,
    AutoGroupsMarkersMissing,
    AutoGroupsParseError,
    AutoGroupsConfigInvalid,
    AutoGroupsAntiParse,
    AutoGroupsEmptyAssignment,
}

impl ErrorKind {
    pub fn code(self) -> &'static str {
        match self {
            Self::InputMissing => "E_INPUT_MISSING",
            Self::TsvHeader => "E_TSV_HEADER",
            Self::TsvParse => "E_TSV_PARSE",
            Self::Path => "E_PATH",
            Self::InvalidArgument => "E_INVALID_ARGUMENT",
            Self::ResourcesIncomplete => "E_RESOURCES_INCOMPLETE",
            Self::GroupsFormat => "E_GROUPS_FORMAT",
            Self::SecretionInputMissing => "E_SECRETION_INPUT_MISSING",
            Self::RegimeMapMissing => "E_REGIME_MAP_MISSING",
            Self::RegimeMapParse => "E_REGIME_MAP_PARSE",
            Self::AutoGroupsMarkersMissing => "E_AUTO_GROUPS_MARKERS_MISSING",
            Self::AutoGroupsParseError => "E_AUTO_GROUPS_PARSE_ERROR",
            Self::AutoGroupsConfigInvalid => "E_AUTO_GROUPS_CONFIG_INVALID",
            Self::AutoGroupsAntiParse => "E_AUTO_GROUPS_ANTI_PARSE",
            Self::AutoGroupsEmptyAssignment => "E_AUTO_GROUPS_EMPTY_ASSIGNMENT",
        }
    }

    pub fn exit_code(self) -> i32 {
        match self {
            Self::InputMissing => 2,
            Self::TsvHeader => 3,
            Self::TsvParse => 4,
            Self::Path => 5,
            Self::InvalidArgument => 6,
            Self::ResourcesIncomplete => 7,
            Self::GroupsFormat => 8,
            Self::SecretionInputMissing => 9,
            Self::RegimeMapMissing => 10,
            Self::RegimeMapParse => 11,
            Self::AutoGroupsMarkersMissing => 12,
            Self::AutoGroupsParseError => 13,
            Self::AutoGroupsConfigInvalid => 14,
            Self::AutoGroupsAntiParse => 15,
            Self::AutoGroupsEmptyAssignment => 16,
        }
    }
}

#[derive(Debug)]
pub struct KiraError {
    pub kind: ErrorKind,
    pub message: String,
}

impl KiraError {
    pub fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
        }
    }
}

impl Display for KiraError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.kind.code(), self.message)
    }
}

impl std::error::Error for KiraError {}

pub type Result<T> = std::result::Result<T, KiraError>;
