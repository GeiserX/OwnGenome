//! What the caller passes on stdin: `{"results": <the results body the page
//! gets>}`. Phase 1 reads the person's sex from it: the sex the reads show
//! (`person.inferred_sex`) when it is male or female, else the declared one.

use pycompat::json::{self, Value};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    /// Neither measured nor declared: no flag that depends on sex is set.
    Unknown,
}

impl Sex {
    pub fn word(self) -> &'static str {
        match self {
            Sex::Male => "male",
            Sex::Female => "female",
            Sex::Unknown => "unknown",
        }
    }

    fn from_word(s: &str) -> Sex {
        match s.to_ascii_lowercase().as_str() {
            "male" => Sex::Male,
            "female" => Sex::Female,
            _ => Sex::Unknown,
        }
    }
}

/// The facts the map needs from stdin.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Facts {
    pub sex: Sex,
}

/// Parses the stdin body. Empty input is no facts; anything else must be a
/// JSON object.
pub fn parse(bytes: &[u8]) -> Result<Facts, String> {
    let text = std::str::from_utf8(bytes).map_err(|_| "stdin is not UTF-8".to_owned())?;
    if text.trim().is_empty() {
        return Ok(Facts { sex: Sex::Unknown });
    }
    let v = json::loads(text).map_err(|e| format!("stdin is not JSON: {e}"))?;
    let Value::Object(top) = v else {
        return Err("stdin is not a JSON object".to_owned());
    };
    let person = match top.get("results") {
        Some(Value::Object(r)) => r.get("person"),
        _ => None,
    };
    let word = |key: &str| match person {
        Some(Value::Object(p)) => match p.get(key) {
            Some(Value::Str(s)) => Sex::from_word(s),
            _ => Sex::Unknown,
        },
        _ => Sex::Unknown,
    };
    let sex = match word("inferred_sex") {
        Sex::Unknown => word("declared_sex"),
        known => known,
    };
    Ok(Facts { sex })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sex_from_reads_then_declared() {
        let p = |s: &str| parse(s.as_bytes()).map(|f| f.sex);
        assert_eq!(p(""), Ok(Sex::Unknown));
        assert_eq!(p("{}"), Ok(Sex::Unknown));
        assert_eq!(
            p(r#"{"results":{"person":{"inferred_sex":"female","declared_sex":"male"}}}"#),
            Ok(Sex::Female)
        );
        assert_eq!(
            p(r#"{"results":{"person":{"inferred_sex":"unknown","declared_sex":"Male"}}}"#),
            Ok(Sex::Male)
        );
        assert_eq!(
            p(r#"{"results":{"person":{"inferred_sex":null,"declared_sex":null}}}"#),
            Ok(Sex::Unknown)
        );
        assert!(p("[1]").is_err());
        assert!(p("{").is_err());
    }
}
